#! /bin/tcsh -f
#
# post-release changes (27 Oct 2026, claude_numbers session; base = 1aho/compare_to_best_runme.com of 5 Jul 2026, the LGD final driver):
#  - this script's own directory goes FIRST on PATH, so helpers kept next to it are used (before, 1aho/ came first and silently
#    overrode them);
#  - uses flip_to_target.csh (UnTangle repo; Arg NH1/NH2 named canonically from chi5) when it sits next to this script, else
#    flip_to_target_runme.com from PATH;
#  - uses swap_compare built next to this script (gcc -O2 -o swap_compare swap_compare.c -lm) when present, else the old path;
#  - logs which phenix.refine, flip script and swap engine it uses.
#  Other absolute paths (1aho data, amber.old/1aho_refine helpers, ccp4) are unchanged.  The paper's LGD runs used phenix 1.21rc1;
#  note that a login .cshrc that sources a phenix_env can override the PATH this script inherits.
#
# changes 29 Oct 2026 (claude_numbers session; JH: "update errant copies ... add set path = ( `dirname $0` $path )",
# "implement fixed-A/B frame at round 0", "implement escalation to larger fragments if single-atom move finder fails"):
#  - set path = ( this script's own directory, absolute, $path ) after every other PATH change, and the per-round comparison
#    calls this script itself (${mypath}/${myname}), not whatever compare_to_best_runme.com comes first on PATH;
#  - phenix_bin=<dir> (default /programs/phenix-1.21rc1-4903/build/bin when it exists; phenix_bin=none to keep PATH as is):
#    strips other phenix entries from PATH and puts this one first (the paper's LGD runs used 1.21rc1);
#  - partition=<name> (default debug): srun partition; the pre-refinement now also runs through $srun, so its result does not
#    depend on the host the driver was started on (29 Oct test: the head node bl831 gives a 0.004 A rms different pre-refined
#    model than every refmac node, and that alone decided whether c-a-c fused Arg18 got trapped);
#  - fixframe=1 (passed by the swap-down loop to every per-round comparison unless keepframe=0): the global A/B frame chosen
#    in round 0 is kept.
#    Round 0 picks it as before (RMSD among confsplit / flip2match / flip2matchBA, then needed-swap count AB vs BA); later rounds
#    may still choose confsplit or flip2match (renaming only) but never flip2matchBA, and AB is forced. Before, a model with about
#    half its atoms in each frame (phenix_ensemble) switched frame after move 1, and move 2 undid move 1;
#  - number_of_needed_swaps.txt is cleared before each comparison (the engine appends; the old AB decision took the minimum over
#    all rounds so far);
#  - pullin=1 (default): when no single-atom A/B relabel is needed but the model has not arrived (a heavy, non-water atom more than
#    arrive_tol = 0.25 A from the target after renaming flips), make a pull-in move: put the worst atom (then larger fragments)
#    at the target position, hold it through the lo and hi jobs, refine as for a swap, compare again. If the same atom is still
#    the worst afterwards, escalate: atom -> side chain (this altloc) -> residue (this altloc) -> residue (both altlocs) ->
#    residues i-1..i+1 -> i-2..i+2 -> every off atom -> stuck. Productive pull-in moves are counted separately from swaps;
#  - the end of the run reports arrived (maxD) or stuck (atoms still off), not just "all done".
#
# changes 30 Oct 2026 (claude_numbers session; JH: "make always_hold apply to every lucky guess move. should be no repeats"):
#  - holdall=1 (default; holdall=0 = the 29 Oct behaviour): every atom a move fixes goes on always_hold.txt, i.e. after each
#    swap move, the swapped atoms that no longer need a swap, and after each pull-in move, the pulled heavy atoms. Held atoms
#    are restrained (reference coordinates, sigma 0.01 A) in every later lo and hi job; the normal refinement of each round
#    stays unrestrained, as before;
#  - a held atom is never chosen as a move again: the move is the top of needed_swaps.txt minus always_hold.txt
#    (needed_new.txt). A held atom that is mislabelled again is re-swapped with the next move without being counted (as the
#    July always_hold did), and logged in silent_reswaps.txt; when only held atoms need a swap, a re-swap round runs on its own
#    (not counted; more than 3 in a row = out of ideas). A pull-in of an atom pulled before is not counted either (re-pull);
#  - the July promotion of struggling atoms to always_hold is off under holdall (an atom whose move never worked would then
#    be re-swapped silently every round instead of escalating);
#  - moves.txt: one line per round (round, kind, counted or not, atoms), so repeats can be checked.
#
#
# check a provided model against the right answer
#
#
set pdbfile = ""
set swapdown = 0

set weight_snap = 0
set hold_first = 1 
set fast = 0
set prerefine = 1
set molprobify = 1
set gt_water = 1
set srun = "srun -p debug"
set debug = 0
set partition = debug
set phenix_bin = /programs/phenix-1.21rc1-4903/build/bin
set fixframe = 0
set keepframe = 1
set pullin = 1
set arrive_tol = 0.25
set holdall = 1

# read the command line to update variables and other settings
foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$key" == "debug") set debug = "1"
end

if(! -e "$pdbfile") then
    set BAD = "pdb file $pdbfile does not exist."
    goto exit
endif

set srun = "srun -p $partition"
if( "$phenix_bin" != "none" && -d "$phenix_bin" ) then
  setenv PATH `echo $PATH | tr ':' '\n' | grep -v phenix | tr '\n' ':'`
  setenv PATH ${phenix_bin}:${PATH}
  unsetenv PHENIX PHENIX_PREFIX PHENIX_VERSION
  rehash
endif

cat << EOF
pdbfile = $pdbfile
weight_snap = $weight_snap
fast = $fast
swapdown = $swapdown
gt_water = $gt_water
debug = $debug
partition = $partition
fixframe = $fixframe
keepframe = $keepframe
pullin = $pullin
arrive_tol = $arrive_tol
holdall = $holdall
EOF
which phenix.refine


set pdir = /home/jamesh/projects/build_challenge/1aho/
set mypath = `dirname $0`
set mypath = `cd $mypath ; pwd`
set myname = `basename $0`
setenv PATH /home/jamesh/projects/amber.old/1aho_refine:${pdir}:${PATH}
# this script's own directory goes first, so the helpers next to it (and this script itself) are the ones used
set path = ( $mypath $path )
echo "driver: ${mypath}/${myname}"
set flipscript = flip_to_target_runme.com
if( -x ${mypath}/flip_to_target.csh ) set flipscript = ${mypath}/flip_to_target.csh
echo "flip script: `which $flipscript`"
rehash

#cp -p ../ground_factory5/nonoise.mtz refme.mtz
cp -p ${pdir}/simplenoise1/noise_borrow0.mtz refme.mtz
cp -p ${pdir}/ground_factory5/fmodelme.pdb ground_truth.pdb

if(! -e gt.pdb && -e ${pdir}/ground_factory5/onechain.pdb) then
  cp ${pdir}/ground_factory5/onechain.pdb gt.pdb
endif
cp -p ${pdir}/refine_best/betteropts_001.pdb good.pdb

cp -p ${pdir}/better_opts.eff opts.eff

# first test convention?
#echo "enforcing"
cat $pdbfile |\
awk '! /^ATOM|^HETAT/{print;next}\
  /HOH/{print substr($0,1,21)"S"substr($0,23);next}\
  {print substr($0,1,21)"A"substr($0,23)}' |\
cat >! provided_onechain.pdb

# canonicalize N-terminal H to refmac naming ("H") to match the refmac-refined
# references (good.pdb / gt.pdb); keeps internal matching consistent and leaves the
# line-201 H1->H duplication dormant.  phenix consumers convert to H1 as needed.
fix_nterm_H.awk -v output=refmac provided_onechain.pdb >! provided_Hfix.pdb
mv provided_Hfix.pdb provided_onechain.pdb

#rmsd provided.pdb gt.pdb | tee rmsd_provided_vs_true.txt
echo "provided vs best"
filter_pdb.awk -v skip=H,water provided_onechain.pdb good.pdb |\
 rmsd >! rmsd_provided_vs_best.txt
head rmsd_provided_vs_best.txt
set warnings = `grep WARNING rmsd_provided_vs_best.txt | wc -l`
set nothing = `grep "no atom pairs" rmsd_provided_vs_best.txt | wc -l`
if ( $nothing ) set warnings = "100%"
if( "$warnings" != "0" ) then
    echo "WARNING: $warnings name changes"
    # analyze name changes
    combine_pdbs_runme.com provided_onechain.pdb good.pdb xor=1 outfile=xor.pdb >! xor.log
    egrep "^ATOM|^HETAT" xor.pdb >! gt_miss.pdb
    combine_pdbs_runme.com good.pdb provided_onechain.pdb xor=1 outfile=xor.pdb >! xor.log
    egrep "^ATOM|^HETAT" xor.pdb >! prov_extra.pdb

    set gt_miss = `cat gt_miss.pdb | wc -l`
    set gt_miss_noalt = `awk 'substr($0,17,1)==" "' provided_onechain.pdb | wc -l`
    set gt_miss_pro = `grep -v HOH gt_miss.pdb | awk '$NF!="H"' | wc -l`
    set gt_miss_H = `awk '$NF=="H"' gt_miss.pdb | wc -l`
    if( $gt_miss ) echo "WARNING: $gt_miss atoms missing from provided"
    if( $gt_miss_noalt ) echo "WARNING: $gt_miss_noalt atoms in provided have no altlocs"
    if( $gt_miss_pro ) echo "WARNING: $gt_miss protein atoms missing from provided"
    if( $gt_miss_H ) echo "WARNING: $gt_miss_H H atoms missing from provided"

    set prov_extra = `cat prov_extra.pdb | wc -l`
    set prov_miss_noalt = `awk 'substr($0,17,1)==" "' prov_extra.pdb | wc -l`
    set prov_extra_pro = `grep -v HOH prov_extra.pdb | awk '$NF!="H"' | wc -l`
    set prov_extra_H = `awk '$NF=="H"' prov_extra.pdb | wc -l`
    if( $prov_extra ) echo "WARNING: $prov_extra extra atoms in provided"
    if( $prov_extra_pro ) echo "WARNING: $prov_extra extra protein atoms in provided"
    if( $prov_extra_H ) echo "WARNING: $prov_extra_H extra H atoms in provided"
endif
echo ""


# split any merged conformers
echo "splitting any merged conformers"
cat provided_onechain.pdb |\
awk '! /^ATOM|^HETAT/{print;next}\
  {conf=substr($0,17,1);occ=substr($0,55,6);\
   X=substr($0,31,8);YZ=substr($0,39,16);\
   pre=substr($0,1,16);typ=substr($0,18,13);post=substr($0,61)}\
  conf==" "{occ=sprintf("%6.2f",occ/2);\
  print pre "A" typ X YZ occ post;\
  X=sprintf("%8.3f",X+0.001);\
  print pre "B" typ X YZ occ post;next}\
  {print}' |\
cat >! confsplit.pdb

# now find rotamer flip/ AB-swap combination that best matches ground truth

$flipscript confsplit.pdb gt.pdb >! ftt.log
mv flipped.pdb flip2match.pdb

swapAB.awk -v swapchain=0 confsplit.pdb >! swap_BA.pdb
$flipscript swap_BA.pdb gt.pdb >! ftt.log
mv flipped.pdb flip2matchBA.pdb

filter_pdb.awk -v skip=H,water confsplit.pdb good.pdb |\
 rmsd >! rmsd_confsplit_vs_best.txt
head rmsd_confsplit_vs_best.txt
set warnings = `grep WARNING rmsd_confsplit_vs_best.txt | wc -l`
set nothing = `grep "no atom pairs" rmsd_confsplit_vs_best.txt | wc -l`
if ( $nothing ) set warnings = "100%"
if( "$warnings" != "0" ) then
    echo "WARNING: $warnings changes after confsplit"
    # analyze name changes
    combine_pdbs_runme.com confsplit.pdb good.pdb xor=1 outfile=xor.pdb >! xor.log
    egrep "^ATOM|^HETAT" xor.pdb >! gt_miss.pdb
    combine_pdbs_runme.com good.pdb confsplit.pdb xor=1 outfile=xor.pdb >! xor.log
    egrep "^ATOM|^HETAT" xor.pdb >! prov_extra.pdb

    set gt_miss = `cat gt_miss.pdb | wc -l`
    set gt_miss_noalt = `awk 'substr($0,17,1)==" "' provided_onechain.pdb | wc -l`
    set gt_miss_pro = `grep -v HOH gt_miss.pdb | awk '$NF!="H"' | wc -l`
    set gt_miss_H = `awk '$NF=="H"' gt_miss.pdb | wc -l`
    if( $gt_miss ) echo "WARNING: $gt_miss atoms missing from confsplit"
    if( $gt_miss_noalt ) echo "WARNING: $gt_miss_noalt atoms in confsplit have no altlocs"
    if( $gt_miss_pro ) echo "WARNING: $gt_miss protein atoms missing from confsplit"
    if( $gt_miss_H ) echo "WARNING: $gt_miss_H H atoms missing from confsplit"

    set prov_extra = `cat prov_extra.pdb | wc -l`
    set prov_miss_noalt = `awk 'substr($0,17,1)==" "' prov_extra.pdb | wc -l`
    set prov_extra_pro = `grep -v HOH prov_extra.pdb | awk '$NF!="H"' | wc -l`
    set prov_extra_H = `awk '$NF=="H"' prov_extra.pdb | wc -l`
    if( $prov_extra ) echo "WARNING: $prov_extra extra atoms in confsplit"
    if( $prov_extra_pro ) echo "WARNING: $prov_extra extra protein atoms in confsplit"
    if( $prov_extra_H ) echo "WARNING: $prov_extra_H extra H atoms in confsplit"
endif

echo ""

echo "rotamer-flipped to target:"
filter_pdb.awk -v skip=H,water flip2match.pdb good.pdb |\
 rmsd >! rmsd_flip2match_vs_best.txt
head rmsd_flip2match_vs_best.txt
echo ""

echo "rotamer-flipped and A-B swapped:"
filter_pdb.awk -v skip=H,water flip2matchBA.pdb good.pdb |\
 rmsd >! rmsd_flip2matchBA_vs_best.txt
head rmsd_flip2matchBA_vs_best.txt
echo ""


egrep RMSD.all rmsd_*_vs_best.txt | sort -k2g | tee sorted.txt | awk -F ":" '{print $2,$1}'
set prefix = `awk -F "_" '{print $2;exit}' sorted.txt | egrep -v "provided|gtHOH"`
if( $fixframe ) then
  # keep the global A/B frame chosen in round 0: renaming flips allowed, no global A<->B swap
  grep -v flip2matchBA sorted.txt >! sorted_fixedframe.txt
  set prefix = `awk -F "_" '{print $2;exit}' sorted_fixedframe.txt | egrep -v "provided|gtHOH"`
  if( "$prefix" == "" ) set prefix = confsplit
  echo "fixed frame: flip2matchBA not allowed"
endif
echo "selecting ${prefix}.pdb"
cp ${prefix}.pdb completeme.pdb
echo ""


# get right atom constellation, drop from provided and add from ground truth if needed
echo "adding any missing atoms"
phenix.reduce completeme.pdb |& awk '/^CRYST|^ATOM|^HETAT/' |\
awk '{print} / H1 /{gsub(" H1 "," H  ");print}' >! completemeH.pdb
egrep -v HOH gt.pdb >! dry_gt.pdb
combine_pdbs_runme.com completemeH.pdb dry_gt.pdb printref=1 saveOcc=1 outfile=dry.pdb >! ftt.log
combine_pdbs_runme.com completemeH.pdb dry_gt.pdb      xor=1 outfile=xor.pdb >> ftt.log
egrep "^ATOM|^HETAT" xor.pdb >! gt_bleed.pdb

set gt_added = `cat gt_bleed.pdb | wc -l`
set gt_added_HOH = `grep HOH gt_bleed.pdb | wc -l`
set gt_added_pro = `grep -v HOH gt_bleed.pdb | awk '$NF!="H"' | wc -l`
set gt_added_H = `awk '$NF=="H"' gt_bleed.pdb | wc -l`

if( $gt_added ) then
    echo "WARNING: $gt_added ground truth atoms added to test model"
else
    echo "no ground truth protein atoms in test model"
endif
if( $gt_added_pro ) echo "WARNING: added $gt_added protein atoms from ground truth that were missing"
if( $gt_added_H ) echo "WARNING: added $gt_added_H H atoms from ground truth that were missing"
if( $gt_added_HOH ) echo "WARNING: added $gt_added_HOH water atoms from ground truth that were missing"


if( $gt_added ) then
    # do something to remove bias
    #combine_pdbs_runme.com confsplit.pdb gt.pdb xor=1 outfile=xor.pdb > /dev/null
    echo "geo-minimizing added ground-truth atoms"
    cat xor.pdb |\
    awk '/^ATOM|^HETAT/{print substr($0,12,5),substr($0,17,1),substr($0,22,1),substr($0,23,6)}' |\
    awk '{print "name",$1,"and altloc",$2,"and chain",$3,"and resid",$4}' |\
    awk '{s = $0}\
        {printf("(%s) or ",s)}' |\
    awk '{$NF="";print}' |\
    awk '{print "selection = \"" $0 "\""}' >! selection.eff
    # phenix.geometry_minimization needs the N-terminal H as H1 or it drifts
    fix_nterm_H.awk -v output=phenix dry.pdb >! dry_Hfix.pdb
    mv dry_Hfix.pdb dry.pdb
    phenix.geometry_minimization dry.pdb \
      selection.eff \
      stop_for_unknowns=false \
      flip_symmetric_amino_acids=False \
      const_shrink_donor_acceptor=0.6 \
      allow_polymer_cross_special_position=True \
        automatic_linking.amino_acid_bond_cutoff=2.5 \
        automatic_linking.inter_residue_bond_cutoff=3.5 \
      output_file_name_prefix=debias >! debias.log
    filter_pdb.awk -v skip=H debias.pdb dry.pdb | rmsd | grep -v Bfac | grep MAXD
    awk '$NF=="H"' debias.pdb dry.pdb | rmsd | grep -v Bfac | grep MAXD
    cp debias.pdb dry.pdb
endif



# now deal with water
filter_pdb.awk -v skip=water dry.pdb | egrep -v "END" >! new.pdb
if( $gt_water ) then
  echo "replacing water with ground truth water"
  filter_pdb.awk -v only=water,atoms good.pdb >> new.pdb
  echo ""
else
  echo "using water from ${prefix}.pdb"
  filter_pdb.awk -v only=water,atoms ${prefix}.pdb >> new.pdb
endif
mv new.pdb refme.pdb

# phenix needs the N-terminal H named H1 (refmac/reduce call it H) or geometry drifts
fix_nterm_H.awk -v output=phenix refme.pdb >! refme_Hfix.pdb
mv refme_Hfix.pdb refme.pdb

if( $prerefine ) then
  if( -e optref_001.pdb && $debug ) then
    echo "DEBUG: re-using optref_001.pdb"
    rmsd optref_001.pdb refme.pdb
  else
    echo "pre-refining with same options as good.pdb"
    set cyc = 10
    if( $fast > 1 ) set cyc = 1
    $srun phenix.refine opts.eff refme.mtz refme.pdb prefix=optref main.number_of_macro_cycles=$cyc >! phenix1.log
  endif
  set optrefpdbfile = optref_001.pdb
else
  echo "using ${prefix}.pdb as-is"
  set optrefpdbfile = refme.pdb
endif

set reffile = good.pdb



swapAB.awk -v swapchain=0 $optrefpdbfile >! BA.pdb
cp $optrefpdbfile AB.pdb

# ===== fast in-memory C reimplementation of the per-atom swap/compare engine =====
# Replaces the original `foreach AB ( BA AB )` loop below (which spawned ~12,000
# flip_to_target/rmsd/awk processes per run). Produces byte-identical
# needed_swaps_*.txt / number_of_needed_swaps.txt and numerically-identical
# sortme_*.txt / min_rmsds_*.txt. Source: claude/swap_compare.c
set swapengine = /home/jamesh/projects/build_challenge/1aho/claude/swap_compare
if( -x ${mypath}/swap_compare ) set swapengine = ${mypath}/swap_compare
echo "swap engine: $swapengine"
if(! -x "$swapengine") then
  echo "ERROR: swap_compare engine not found/executable at $swapengine"
  echo "build it with: gcc -O2 -o $swapengine ${swapengine}.c -lm"
  exit 9
endif
rm -f number_of_needed_swaps.txt
$swapengine good.pdb BA AB

# ----- original awk engine kept inline for reference (disabled) -----
if ( 0 ) then
foreach AB ( BA AB )
  #set AB = AB

  rm -f min_rmsds.txt
  rm -f all_rmsds.txt
  rm -f bconfs.txt

  cat ${AB}.pdb |\
  awk '$NF!="H" && /^ATOM|^HETAT/{\
    atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
    gsub(" ","",atom);print atom,chain resnum}' |\
   sort -u | sort -k2.1,2.1 -k2.2g |\
  cat >! atoms.txt
  set atoms = `awk '{print $1"_"$2}' atoms.txt`

  rm -f sortme.txt
  foreach atom ( $atoms )

    set cr = `echo $atom | awk -F "_" '{print $2}'`

    cat ${AB}.pdb |\
    awk -v cr=$cr '$NF=="H"{next} ! /^ATOM|^HETAT/{next}\
      {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
      gsub(" ","",atom);this=chain resnum}\
      this==cr{print}' >! bres.pdb 
    cat $reffile |\
    awk -v cr=$cr '$NF=="H"{next} ! /^ATOM|^HETAT/{next}\
      {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
      gsub(" ","",atom);this=chain resnum}\
      this==cr{print}' >! gres.pdb 

    cat bres.pdb |\
    awk '! /^ATOM|^HETAT/{print;next}\
      {c=substr($0,17,1);pre=substr($0,1,16);post=substr($0,18)}\
      c==" "{print;next}\
      c=="A"{print pre"B"post}\
      c=="B"{print pre"A"post}' |\
    cat >! sres.pdb

    $flipscript bres.pdb gres.pdb >! bflip.log
    mv flipped.pdb bres.pdb
    $flipscript sres.pdb gres.pdb >! sflip.log
    mv flipped.pdb sres.pdb
    set bflip = `egrep MAXD.all bflip.log | awk '$2==0{exit} {print}' | wc -l`
    set sflip = `egrep MAXD.all sflip.log | awk '$2==0{exit} {print}' | wc -l`

    cat bres.pdb |\
    awk -v a=$atom '$NF=="H"{next} ! /^ATOM|^HETAT/{next}\
      {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
      gsub(" ","",atom);this=atom"_"chain resnum}\
      this==a{print}' >! batm.pdb 
    cat gres.pdb |\
    awk -v a=$atom '$NF=="H"{next} ! /^ATOM|^HETAT/{next}\
      {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
      gsub(" ","",atom);this=atom"_"chain resnum}\
      this==a{print}' >! gatm.pdb 
    cat sres.pdb |\
    awk -v a=$atom '$NF=="H"{next} ! /^ATOM|^HETAT/{next}\
      {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
      gsub(" ","",atom);this=atom"_"chain resnum}\
      this==a{print}' >! satm.pdb 

    set nstat = `awk '$NF!="H"' gatm.pdb batm.pdb | rmsd -v xlog=1 | awk '{print $2}'`
    set sstat = `awk '$NF!="H"' gatm.pdb satm.pdb | rmsd -v xlog=1 | awk '{print $2}'`
    set dstat = `awk '$NF!="H"' batm.pdb satm.pdb | rmsd -v xlog=1 | awk '{print $2}'`

    echo "$atom $sstat $nstat $dstat  $sflip $bflip" | tee -a sortme.txt > /dev/null
    if( $debug ) tail -n 1 sortme.txt
  end
  # atomcode good-swap good-bad delta-swap  sflip bflip
  #    1         2        3         4         5      6

  cp sortme.txt sortme_${AB}.txt

  cat sortme.txt |\
  awk '{print $1,$2,"B",$5;print $1,$3,"A",$6}' |\
  sort -k2g |\
  awk '! seen[$1]{print;++seen[$1]}' |\
  cat >! min_rmsds.txt
  cp min_rmsds.txt min_rmsds_${AB}.txt
  # atomcode min_rmsd conf flip 
  if( $debug ) cat min_rmsds.txt

  grep -v _S sortme_${AB}.txt |\
  awk '$2<$3{print ($3-$2)**2+$4**2+$3**2,$0}' |\
  sort -gr |\
  awk '{print "SWAP",$2,$1,$3,$4,$5}' >! needed_swaps_${AB}.txt
  # SWAP atomcode score good-swap good-bad delta-swap

  set n = `wc -l needed_swaps_${AB}.txt | awk '{print $1}'`

  echo "$n $AB" | tee -a number_of_needed_swaps.txt

end
endif
# end of AB loop (original awk engine; replaced by $swapengine above)

set AB = `sort -g number_of_needed_swaps.txt | awk '{print $2;exit}'`
if( $fixframe ) then
  cat number_of_needed_swaps.txt
  set AB = AB
  echo "fixed frame: AB kept"
endif
echo "decided on AB = $AB"
cp sortme_${AB}.txt sortme.txt
cp needed_swaps_${AB}.txt needed_swaps.txt

awk '$NF==1' min_rmsds.txt |\
awk '{print $1}' | awk -F "_" '{print $2}' |\
sort -u | sort -k1.2g >! unique_minrmsds.txt
if( $debug ) cat unique_minrmsds.txt

cp ${AB}.pdb bad.pdb


awk '{print $2,sprintf("%.1f",$3)}' needed_swaps.txt |\
  awk -F "[_ ]" '! worst[$2]{worst[$2]=$3} {atoms[$2]=atoms[$2]" "$1}\
  END{for(r in atoms)print r,worst[r],atoms[r]}' |\
 sort -k2gr >! swaps_byres.txt

if( $molprobify ) then
set pfx = bad
set good = good

echo "molprobify both..."
molprobify_runme.com ${good}.pdb keepgeo >&! molprobify_best.log &
molprobify_runme.com ${pfx}.pdb keepgeo >&! molprobify_${pfx}.log
$flipscript ${pfx}.pdb good.pdb 
wait
tail -n 1 molprobify_best.log | awk '{print $(NF-2) ":","wE:",$2,"R:",$3,"Rfree:",$4}'
tail -n 1 molprobify_${pfx}.log | awk '{print $(NF-2) ":","wE:",$2,"R:",$3,"Rfree:",$4}'

foreach prefix ( ${good} ${pfx} )
  cat ${prefix}_fullgeo.txt |\
  awk -F "|" '{n=split($2,w,/[ \t]+/);\
   #print;print w[3],w[9],w[15],w[21];\
   w[3]=w[9]=w[15]=w[21]="";\
   s="";for(i=1;i<=n;++i)s=s" "w[i];\
   print $1"|",s}' |\
  awk -F "|" '! seen[$2]{print;++seen[$2]}' |\
  cat >! ${prefix}_worst.txt
end


awk '{print $0 "|GOOD"}' ${good}_worst.txt |\
 cat - ${pfx}_worst.txt | awk '{$1=$1"|";print}' |\
  awk -F "|" '{id=$1" |"$3;gsub(/ $|^ /,"",id);++count[id]}\
  $NF=="GOOD"{good[id]=$2+0;next}\
 $2+0>0{print $2-good[id],"sigma:",sqrt(good[id]*(good[id]>0)),"->",sqrt($2+0),"for",id}' |\
tee diffs.txt | wc -l

echo "dwE          good -> bad "
sort -gr diffs.txt | head | awk '{$2="";print}' | justify.awk 
sort -gr diffs.txt | tail | awk '{$2="";print}' | justify.awk 

awk '$NF!="H" && ! /HOH/' flipped.pdb good.pdb | rmsd -v debug=1 |\
 sort -k1.25g | grep moved |\
 awk '{gsub("moved","maxD:");\
  print substr($0,1,5),substr($0,6,1),substr($0,7,28)}' |\
 awk '{print $6,$7,$1,$2,$3,$4,$5}' |\
  tail -n 1 >! maxd.txt
set maxd = `tail -n 1 maxd.txt`

set worstdelta = `sort -g diffs.txt | tail -n 1 | awk '{gsub("\\|"," ");print}'`
set bestdelta = `sort -g diffs.txt | head -n 1 | awk '{gsub("\\|"," ");print}'`

#echo "dwE  good -> bad "
#sort -gr diffs.txt | head | awk '{$2="";print}' | justify.awk 
#sort -gr diffs.txt | tail | awk '{$2="";print}' | justify.awk 
endif
# end of molprobify


echo "swaps by residue:"
cat swaps_byres.txt

if( ! $swapdown ) exit

# 
swapdown:

echo "entering swap-down phase..."
set productive_steps = 0

# save from being overwritten
mkdir -p before_swaps
cp -p *.pdb *.txt before_swaps
cp -p *.mtz *.log before_swaps
cp -p *.geo       before_swaps

cp $pdbfile provided.pdb
cp $optrefpdbfile provided_optref.pdb
molprobify_runme.com provided.pdb >! molprobify_provided.log &
molprobify_runme.com provided_optref.pdb >! molprobify_provided_optref.log &


if(! -e Hbonds.txt) then
  gemmi rmsz --cutoff=-1 good.pdb |\
   tee gemmi_rmsz.txt |\
   awk '$3=="bond" && $4 ~/-H/{print $4}' |\
  awk -F "[-:]" '{print "BOND",$1,$2}' |\
  cat >! Hbonds.txt

endif

set swaps = 0

if( $gt_water ) then
  echo "replacing water with ground truth water"
  filter_pdb.awk -v skip=water bad.pdb | egrep -v "END" >! swapme.pdb
  filter_pdb.awk -v only=water good.pdb >> swapme.pdb
else
   cp bad.pdb swapme.pdb
endif


cp molprobify_bad.log molprobify_oaat_${swaps}.log

set prev = ""
echo -n "" >! needed_swaps_${prev}.txt
echo -n "" >! needed_new_${prev}.txt

# atoms whose swap required holdme.eff to secure (one-shot list);
# any of those that reappear in needed_swaps.txt get promoted to always_hold
# and are then swap-corrected + held in every subsequent refinement
# (30 Oct: under holdall=1 every atom a move fixes goes on always_hold instead, and there is no promotion)
if(! -e needs_holding.txt) then
   echo -n "" >! needs_holding.txt
endif
if(! -e always_hold.txt) then
   echo -n "" >! always_hold.txt
endif
if(! -e moves.txt) then
   echo -n "" >! moves.txt
endif

set struggling = 0
set movetype = one_atom
set prev_movetype = $movetype
set pullin_steps = 0
set pull_level = 0
set pull_prev = ""
set stuck = 0
set swaps_done = 0
set reswap_only = 0
set silent_reswaps = 0
set repulls = 0

movephase:
while ( -s needed_swaps.txt )
@ swaps_done = ( $swaps_done + 1 )

cp needed_swaps.txt needed_swaps_${swaps}.txt

# the move is chosen among the atoms not held yet (holdall); held atoms that need a swap again are re-swapped below, uncounted
if( $holdall ) then
   awk 'FILENAME==ARGV[1]{h[$1]=1;next} ! h[$2]' always_hold.txt needed_swaps.txt >! needed_new.txt
else
   cp needed_swaps.txt needed_new.txt
endif
cp needed_new.txt needed_new_${swaps}.txt
set counted = 0

if( ! -s needed_new.txt ) then
   # only held atoms need a swap: re-swap them in a round of their own, not counted
   @ reswap_only = ( $reswap_only + 1 )
   echo "only held atoms need a swap: re-swap round $reswap_only (not counted):" `awk '{print $2}' needed_swaps.txt`
   if( $reswap_only > 3 ) then
      echo "held atoms keep reverting: out of ideas"
      break
   endif
   echo -n "" >! this_swap.txt
   set struggling = 0
   set movetype = reswap
   set prev_movetype = reswap
else
set reswap_only = 0
head -n 1 needed_new.txt | tee this_swap.txt
set atomid = `awk '{print $2;exit}' needed_new.txt`
set atom = `echo $atomid | awk -F "_" '{print $1}'`
set res = `echo $atomid | awk -F "_" '{print $2}'`
set resnum = `echo $res | awk '{print substr($1,2)}'`
set chain = `echo $res | awk '{print substr($1,1,1)}'`
set rigid = `echo $atom | awk '{print ( $1 ~ /[0-9]$/ || /CZ/ || /OH/ )}'`
set partner = `awk '/^SSBOND/{print $4 $5,$7 $8}' good.pdb |  awk '{print ;print $2,$1}' | egrep "^$res" | awk '{print $2}' | sort -u`

set this_atom = `awk '{print $2;exit}' needed_new.txt | sort`
set prev_atom = `awk '{print $2;exit}' needed_new_${prev}.txt | sort`
set prev_movetype = $movetype
if( "$this_atom" == "$prev_atom" ) then
   echo "swap $swaps did not work for $prev_movetype "
   if( $struggling && "$prev_movetype" == "one_atom" ) then
      # rigid atoms try the rigid side chain first; non-rigid atoms (e.g. a
      # disulfide SG, which never matches the rigid test) have no rigid-side
      # step, so go straight to whole-residue. Either way the escalation keeps
      # climbing toward SSbond/all instead of looping forever on one_atom.
      if( $rigid ) then
         echo "swapping rigid side chain this time "
         set movetype = rigid_side
      else
         echo "swapping whole residue this time (non-rigid atom)"
         set movetype = whole_res
      endif
   endif
   if( $struggling && "$prev_movetype" == "rigid_side" ) then
      echo "swapping whole residue this time "
      set movetype = whole_res
   endif
   if( $struggling && "$prev_movetype" == "whole_res" ) then
      echo "swapping both SGs in disulfide"
      set movetype = SSbond
   endif
   if( $struggling && "$prev_movetype" == "SSbond" ) then
      echo "swapping top pair of needed swaps"
      set movetype = top2
   endif
   if( $struggling && "$prev_movetype" == "top2" ) then
      echo "swapping top 3 tripled of needed swaps"
      set movetype = top3
   endif
   if( $struggling && "$prev_movetype" == "top3" ) then
      echo "swapping everything at once"
      set movetype = "all"
   endif
   if( $struggling && "$prev_movetype" == "all" ) then
      echo "out of ideas"
      break
   endif
   echo "struggling with $prev_movetype"
   set struggling = 1
else
   echo "not struggling"
   set struggling = 0
   @ productive_steps = ( $productive_steps + 1 )
   set counted = 1
   set movetype = one_atom
   set prev_movetype = one_atom
endif

if( "$movetype" == "rigid_side" ) then
   echo "swapping rigid side chain of $res"
   cat needed_new.txt |\
   awk 'NR==1{split($2,w,"_");res=w[2]}\
     $2~/[1-9]_/ && $2~ "_"res"$"{print}' | tee this_swap.txt
endif
if( "$movetype" == "whole_res" ) then
    echo "swapping whole residue $res"
    egrep "_${res} " needed_new.txt | tee this_swap.txt
endif
if( "$movetype" == "SSbond" ) then
    echo "swapping whole disulfide: $res - $partner"
    egrep "_${res} |_${partner} " needed_new.txt | tee this_swap.txt
endif
if( "$movetype" == "top2" ) then
    echo "swapping top 2 entries"
    head -n 2 needed_new.txt | tee this_swap.txt
endif
if( "$movetype" == "top3" ) then
    echo "swapping top 3 entries"
    head -n 3 needed_new.txt | tee this_swap.txt
endif
if( "$movetype" == "all" ) then
    echo "swapping all entries"
    cat needed_new.txt | tee this_swap.txt
endif
endif
awk '/^SWAP/{print $2}' this_swap.txt >! _move_atoms.txt
echo -n "" >! _reswap_now.txt

# augment this_swap.txt with any always_hold atoms that are currently mis-labelled
# (i.e. they appear in needed_swaps.txt) so they get swap-corrected this iteration
# (30 Oct: FILENAME test instead of NR==FNR, which misreads an empty first file, as this_swap.txt is in a re-swap round)
if( -s always_hold.txt ) then
   awk '/^SWAP/{print $2}' this_swap.txt >! _curr_swap_set.txt
   awk 'FILENAME==ARGV[1]{seen[$1]=1;next} ! seen[$1]' _curr_swap_set.txt always_hold.txt >! _ah_to_add.txt
   if( -s _ah_to_add.txt ) then
      awk 'FILENAME==ARGV[1]{want[$1]=1;next} /^SWAP/ && want[$2]' _ah_to_add.txt needed_swaps.txt >! _reswap_now.txt
      cat _reswap_now.txt >> this_swap.txt
      awk '{print "always_hold re-swap:",$0}' _reswap_now.txt
   endif
   rm -f _curr_swap_set.txt _ah_to_add.txt
endif

cat Hbonds.txt this_swap.txt swapme.pdb |\
awk 'BEGIN{other["A"]="B";other["B"]="A"}\
  /^BOND/{hatoms[$2]=hatoms[$2]" "$3;next}\
  /^SWAP/{++swapconf[$2];\
    split($2,w,"_");atm=w[1];rid=w[2];split(hatoms[atm],a);\
    for(i in a){++swapconf[a[i]"_"rid];};next}\
  ! /^ATOM|^HETAT|^ANISO/{next}\
  {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
   gsub(" ","",atom);conf=substr($0,17,1);\
   pre=substr($0,1,16);post=substr($0,18);\
   atm=atom"_"chain resnum}\
  ! swapconf[atm]{next}\
  {print pre other[conf] post}' |\
cat >! swapped.pdb


cat << EOF >! holdme.eff    
refinement {
  main {
    nqh_flips=False
  }
  modify_start_model {
    modify {
      flip_symmetric_amino_acids = False
    }   
  }
  pdb_interpretation {
    flip_symmetric_amino_acids = False
    reference_coordinate_restraints {
      enabled = True
      exclude_outliers = False
EOF
( awk '/^SWAP/{print $2}' this_swap.txt ; cat always_hold.txt ) | sort -u |\
awk 'NF>0{split($1,w,"_");\
     atom=w[1];chain=substr(w[2],1,1);resnum=substr(w[2],2);\
     t="(name " atom " and resseq " resnum " and chain " chain ")";\
     sel=(sel=="")? t : sel " or " t}\
     END{if(sel!="") print "      selection = \"" sel "\""}' |\
cat >> holdme.eff
cat << EOF >> holdme.eff
      sigma = 0.01
      limit = 10
    }
  }
}
EOF

echo -n "" >! cheat.pdb 
set cheating = `awk '{print ( $4>0.2 );exit}' this_swap.txt`
if( $cheating ) then
   echo "cheating a bit because swap-to location >0.2 A away from right place"
   combine_pdbs_runme.com swapped.pdb good.pdb saveXYZ outfile=cheat.pdb > /dev/null
endif
combine_pdbs_runme.com swapped.pdb cheat.pdb swapme.pdb keepref=1 outfile=refme.pdb > /dev/null

# phenix needs the N-terminal H named H1 (refmac/reduce call it H) or geometry drifts
fix_nterm_H.awk -v output=phenix refme.pdb >! refme_Hfix.pdb
mv refme_Hfix.pdb refme.pdb

set prev = "$swaps"
set swaps = `echo $swaps | awk '{printf("%03d",$1+1)}'`
# one line per round: round, move kind, counted (1) or not, the move's atoms, then any held atoms re-swapped with it
set moved = `cat _move_atoms.txt`
set reswapped = `awk '{print $2}' _reswap_now.txt`
echo "$swaps swap $movetype counted=$counted move: $moved reswap: $reswapped" >> moves.txt
if( -s _reswap_now.txt ) then
   awk -v r=$swaps '{print r, $2}' _reswap_now.txt >> silent_reswaps.txt
   set nre = `cat _reswap_now.txt | wc -l`
   @ silent_reswaps = ( $silent_reswaps + $nre )
endif


set hold = ""
if( $struggling ) echo "struggling..."
if( $struggling || $cheating || -s always_hold.txt ) then
    echo "holding swapped atoms in place"
    set hold = holdme.eff
    cp holdme.eff holdme_${swaps}.eff
endif
if( $struggling || $cheating ) then
    # record atoms whose swap needed hold this iteration
    awk '/^SWAP/{print $2}' this_swap.txt >> needs_holding.txt
    sort -u needs_holding.txt -o needs_holding.txt
endif

set cyc = ""
if( $fast ) set cyc = "main.number_of_macro_cycles=1"

if ( $weight_snap ) then
  echo "refine wxc_scale=0.05"
  set thishold = "$hold"
  if( $hold_first ) set thishold = holdme.eff
  $srun phenix.refine refme.mtz opts.eff $thishold $cyc \
    automatic_linking.amino_acid_bond_cutoff=2.5 \
    automatic_linking.inter_residue_bond_cutoff=3.5 \
    refme.pdb prefix=lo wxc_scale=0.05 >&! low_${swaps}.log 
  cp lo_001.pdb refme.pdb
  echo "refine wxc_scale=5"
  $srun phenix.refine refme.mtz opts.eff $hold $cyc refme.pdb prefix=hi wxc_scale=5 >&! high_${swaps}.log 
  cp hi_001.pdb refme.pdb
endif
echo "refine normal"
$srun phenix.refine refme.mtz opts.eff $cyc \
    automatic_linking.amino_acid_bond_cutoff=2.5 \
    automatic_linking.inter_residue_bond_cutoff=3.5 \
 refme.pdb prefix=oaat serial=$swaps >&! oaat_${swaps}.log 

cp oaat_${swaps}.pdb swapme.pdb

# save time in comparison, since this refinement has already been done
cp oaat_${swaps}.pdb optref_001.pdb

echo "comparing to best again $swaps ..."
${mypath}/${myname} swapme.pdb prerefine=0 gt_water=$gt_water fixframe=$keepframe partition=$partition phenix_bin=$phenix_bin molprobify=$molprobify >&! compare_to_best_${swaps}.log
# updates lots of files

# holdall: every atom this round swapped that no longer needs a swap is held from now on
if( $holdall ) then
   awk 'FILENAME==ARGV[1]{need[$2]=1;next} /^SWAP/ && ! need[$2]{print $2}' needed_swaps.txt this_swap.txt | sort -u >! _fixed.txt
   if( -s _fixed.txt ) then
      echo "holding from now on:" `cat _fixed.txt`
      cat _fixed.txt always_hold.txt | sort -u >! _ah_new.txt
      mv _ah_new.txt always_hold.txt
   endif
   rm -f _fixed.txt
endif

# any atom that previously needed hold and now reappears in needed_swaps.txt
# is promoted to always_hold (permanent swap-correction + reference restraint); not under holdall
if( ! $holdall && -s needs_holding.txt && -s needed_swaps.txt ) then
   awk 'NR==FNR{want[$1]=1;next} /^SWAP/ && want[$2]{print $2}' \
       needs_holding.txt needed_swaps.txt | sort -u >! _promote.txt
   # only count atoms not already in always_hold
   if( -s always_hold.txt ) then
      awk 'NR==FNR{have[$1]=1;next} ! have[$1]' always_hold.txt _promote.txt >! _new_perm.txt
   else
      cp _promote.txt _new_perm.txt
   endif
   if( -s _new_perm.txt ) then
      echo "promoting to always_hold:" `cat _new_perm.txt | tr '\n' ' '`
      cat _new_perm.txt >> always_hold.txt
      sort -u always_hold.txt -o always_hold.txt
   endif
   rm -f _promote.txt _new_perm.txt
endif

cp molprobify_bad.log molprobify_oaat_${swaps}.log

end

# ===== arrival test and pull-in moves (29 Oct 2026) =====
# Reached when no single-atom A/B relabel is needed (or the swap escalation ran out of ideas).
if( -s needed_swaps.txt ) goto finish

# has the model arrived? heavy, non-water atoms vs the target, after renaming flips (X-ray-equivalent names)
$flipscript swapme.pdb good.pdb >&! arrive_flip.log
cp flipped.pdb arrive_model.pdb
# symmetric names count as equivalent (as label_check_v2): the distance is to the nearer of the atom's own and its partner's
# target position (Asp OD1/OD2, Glu OE1/OE2, Phe/Tyr CD1/CD2 CE1/CE2, Arg NH1/NH2, Val CG1/CG2, Leu CD1/CD2, His ring flip,
# Asn OD1/ND2, Gln OE1/NE2, C-terminal O/OXT)
awk -v tol=$arrive_tol 'BEGIN{n=split("ASP OD1 OD2 GLU OE1 OE2 PHE CD1 CD2 PHE CE1 CE2 TYR CD1 CD2 TYR CE1 CE2 ARG NH1 NH2 VAL CG1 CG2 LEU CD1 CD2 HIS ND1 CD2 HIS CE1 NE2 ASN OD1 ND2 GLN OE1 NE2",w);\
    for(i=1;i<=n;i+=3){pt[w[i]" "w[i+1]]=w[i+2];pt[w[i]" "w[i+2]]=w[i+1]}}\
  ! /^ATOM|^HETAT/{next} /HOH/{next}\
  {name=substr($0,13,4);gsub(" ","",name);el=substr($0,77,2);gsub(" ","",el);rn=substr($0,18,3)}\
  el=="H" || (el=="" && name ~ /^[0-9]*H/){next}\
  {tail="_"substr($0,17,1)"_"substr($0,22,1)"_"(substr($0,23,4)+0);key=name tail;\
   x=substr($0,31,8)+0;y=substr($0,39,8)+0;z=substr($0,47,8)+0}\
  FNR==NR{gx[key]=x;gy[key]=y;gz[key]=z;next}\
  (key in gx){d=sqrt((x-gx[key])^2+(y-gy[key])^2+(z-gz[key])^2);\
   p=pt[rn" "name]; if(p=="" && name=="O")p="OXT"; if(p=="" && name=="OXT")p="O";\
   if(p!="" && (p tail in gx)){d2=sqrt((x-gx[p tail])^2+(y-gy[p tail])^2+(z-gz[p tail])^2);if(d2<d)d=d2}\
   ++n;if(d>max)max=d;if(d>tol)printf("%.3f %s %s\n",d,key,rn)}\
  END{printf("%d %.3f\n",n,max) > "arrive_stats.txt"}' good.pdb arrive_model.pdb |\
  sort -gr >! off_atoms.txt
set arrive_maxd = `awk '{print $2}' arrive_stats.txt`
set nheavy = `awk '{print $1}' arrive_stats.txt`
set noff = `cat off_atoms.txt | wc -l`
echo "arrival test: maxD $arrive_maxd A over $nheavy heavy atoms; $noff more than $arrive_tol A from the target"
if( $noff == 0 ) goto finish
head -n 10 off_atoms.txt
if( ! $pullin ) goto finish

# pull-in move: put the worst atom (or a larger fragment around it) at the target position
set target = `awk '{print $2;exit}' off_atoms.txt`
if(! -e pull_history.txt) touch pull_history.txt
# an atom pulled before starts one level above the highest level already tried on it, so no atom can cycle forever
set triedlevel = `awk -v t=$target '$1==t && $2>m{m=$2} END{print m+0}' pull_history.txt`
set pcounted = 0
if( "$target" == "$pull_prev" ) then
   @ pull_level = ( $pull_level + 1 )
   echo "pull-in of $target did not stick; escalating to level $pull_level"
else
   @ pull_level = ( $triedlevel + 1 )
   set pcounted = 1
   if( $holdall && $triedlevel > 0 ) set pcounted = 0
   if( $pcounted ) then
      set reswap_only = 0
      @ pullin_steps = ( $pullin_steps + 1 )
      if( $triedlevel > 0 ) echo "$target was pulled before (up to level $triedlevel); starting at level $pull_level"
   else
      @ repulls = ( $repulls + 1 )
      echo "$target was pulled before (up to level $triedlevel): re-pull at level $pull_level, not counted (holdall)"
   endif
endif
set pull_prev = "$target"

pullselect:
if( $pull_level > 7 ) then
   echo "out of ideas for $target"
   set stuck = 1
   goto finish
endif
set levelname = ( atom sidechain residue residue_both residues3 residues5 all_off )
echo "pull-in move: $target, level $pull_level ($levelname[$pull_level])"
# Hbonds.txt lists BOND <heavy> <H> by atom name; off_atoms.txt lists every off heavy atom (level 7)
cat Hbonds.txt off_atoms.txt arrive_model.pdb |\
awk -v t=$target -v L=$pull_level 'BEGIN{split(t,w,"_");tn=w[1];ta=w[2];tc=w[3];tr=w[4]+0;\
    n=split("ASP OD1 OD2 GLU OE1 OE2 PHE CD1 CD2 PHE CE1 CE2 TYR CD1 CD2 TYR CE1 CE2 ARG NH1 NH2 VAL CG1 CG2 LEU CD1 CD2 HIS ND1 CD2 HIS CE1 NE2 ASN OD1 ND2 GLN OE1 NE2",q);\
    for(i=1;i<=n;i+=3){pt[q[i]" "q[i+1]]=q[i+2];pt[q[i]" "q[i+2]]=q[i+1]}\
    bb["N"]=bb["CA"]=bb["C"]=bb["O"]=bb["OXT"]=bb["H"]=bb["H1"]=bb["H2"]=bb["H3"]=bb["HA"]=bb["HA2"]=bb["HA3"]=1}\
  /^BOND/{par[$3]=par[$3]" "$2;next}\
  NF==3 && $1+0>0 && $2 ~ /_/{off[$2]=1;next}\
  ! /^ATOM|^HETAT/{next} /HOH/{next}\
  {name=substr($0,13,4);gsub(" ","",name);a=substr($0,17,1);c=substr($0,22,1);r=substr($0,23,4)+0;\
   el=substr($0,77,2);gsub(" ","",el);isH=(el=="H");\
   heavy=name; if(isH){split(par[name],p);heavy="";for(i in p)if(p[i]!="")heavy=heavy" "p[i]}}\
  {sel=0}\
  L==1 && c==tc && r==tr && a==ta{ tp=pt[substr($0,18,3)" "tn]; if(tn=="O")tp="OXT"; if(tn=="OXT")tp="O";\
     if(name==tn || (tp!="" && name==tp)) sel=1;\
     if(isH && (index(" "heavy" "," "tn" ") || (tp!="" && index(" "heavy" "," "tp" ")))) sel=1 }\
  L==2 && c==tc && r==tr && a==ta && ! bb[name]{sel=1}\
  L==3 && c==tc && r==tr && a==ta{sel=1}\
  L==4 && c==tc && r==tr{sel=1}\
  L==5 && c==tc && r>=tr-1 && r<=tr+1{sel=1}\
  L==6 && c==tc && r>=tr-2 && r<=tr+2{sel=1}\
  L==7 && ! isH && off[name"_"a"_"c"_"r]{sel=1}\
  L==7 && isH{n=split(heavy,h);for(i=1;i<=n;++i)if(off[h[i]"_"a"_"c"_"r])sel=1}\
  sel{print}' >! pulled.pdb
if( ! -s pulled.pdb ) then
   echo "nothing to pull at level $pull_level (e.g. Gly side chain); next level"
   @ pull_level = ( $pull_level + 1 )
   goto pullselect
endif
echo "`grep -c . pulled.pdb` atoms pulled in"
echo "$target $pull_level" >> pull_history.txt

cp needed_swaps.txt needed_swaps_${swaps}.txt
combine_pdbs_runme.com pulled.pdb good.pdb saveXYZ outfile=cheat.pdb > /dev/null
combine_pdbs_runme.com cheat.pdb arrive_model.pdb keepref=1 outfile=refme.pdb > /dev/null
fix_nterm_H.awk -v output=phenix refme.pdb >! refme_Hfix.pdb
mv refme_Hfix.pdb refme.pdb

# hold the pulled heavy atoms (and any always_hold atoms) through the lo and hi jobs
cat << EOF >! holdme.eff    
refinement {
  main {
    nqh_flips=False
  }
  modify_start_model {
    modify {
      flip_symmetric_amino_acids = False
    }   
  }
  pdb_interpretation {
    flip_symmetric_amino_acids = False
    reference_coordinate_restraints {
      enabled = True
      exclude_outliers = False
EOF
( awk '{el=substr($0,77,2);gsub(" ","",el)} el!="H"{name=substr($0,13,4);gsub(" ","",name);\
   print "(name " name " and resseq " substr($0,23,4)+0 " and chain " substr($0,22,1) " and altloc " substr($0,17,1) ")"}' pulled.pdb ;\
  awk 'NF>0{split($1,w,"_");print "(name " w[1] " and resseq " substr(w[2],2) " and chain " substr(w[2],1,1) ")"}' always_hold.txt ) |\
awk 'NF>0{sel=(sel=="")? $0 : sel " or " $0} END{if(sel!="") print "      selection = \"" sel "\""}' >> holdme.eff
cat << EOF >> holdme.eff
      sigma = 0.01
      limit = 10
    }
  }
}
EOF

set prev = "$swaps"
set swaps = `echo $swaps | awk '{printf("%03d",$1+1)}'`
cp holdme.eff holdme_${swaps}.eff
( echo "PULLIN $target level $pull_level $levelname[$pull_level]" ; cat pulled.pdb ) >! pullin_${swaps}.txt
echo "$swaps pullin level$pull_level counted=$pcounted move: $target" >> moves.txt

set cyc = ""
if( $fast ) set cyc = "main.number_of_macro_cycles=1"
if ( $weight_snap ) then
  echo "refine wxc_scale=0.05"
  $srun phenix.refine refme.mtz opts.eff holdme.eff $cyc \
    automatic_linking.amino_acid_bond_cutoff=2.5 \
    automatic_linking.inter_residue_bond_cutoff=3.5 \
    refme.pdb prefix=lo wxc_scale=0.05 >&! low_${swaps}.log 
  cp lo_001.pdb refme.pdb
  echo "refine wxc_scale=5"
  $srun phenix.refine refme.mtz opts.eff holdme.eff $cyc refme.pdb prefix=hi wxc_scale=5 >&! high_${swaps}.log 
  cp hi_001.pdb refme.pdb
endif
echo "refine normal"
$srun phenix.refine refme.mtz opts.eff $cyc \
    automatic_linking.amino_acid_bond_cutoff=2.5 \
    automatic_linking.inter_residue_bond_cutoff=3.5 \
 refme.pdb prefix=oaat serial=$swaps >&! oaat_${swaps}.log 

cp oaat_${swaps}.pdb swapme.pdb
cp oaat_${swaps}.pdb optref_001.pdb

echo "comparing to best again $swaps (after pull-in) ..."
${mypath}/${myname} swapme.pdb prerefine=0 gt_water=$gt_water fixframe=$keepframe partition=$partition phenix_bin=$phenix_bin molprobify=$molprobify >&! compare_to_best_${swaps}.log
cp molprobify_bad.log molprobify_oaat_${swaps}.log
# holdall: the pulled heavy atoms (both altlocs of each name, as always_hold selects) are held from now on
if( $holdall ) then
   awk '{el=substr($0,77,2);gsub(" ","",el)} el!="H"{name=substr($0,13,4);gsub(" ","",name);print name "_" substr($0,22,1) (substr($0,23,4)+0)}' pulled.pdb |\
     sort -u >! _fixed.txt
   echo "holding from now on:" `cat _fixed.txt`
   cat _fixed.txt always_hold.txt | sort -u >! _ah_new.txt
   mv _ah_new.txt always_hold.txt
   rm -f _fixed.txt
endif
goto movephase

finish:
set test = `cat needed_swaps.txt | wc -l`
@ total_steps = ( $productive_steps + $pullin_steps )
echo "$test swaps remain"
echo "$productive_steps productive moves in total"
echo "$pullin_steps productive pull-in moves"
echo "$total_steps productive moves including pull-ins"
if( $holdall ) then
  echo "holdall: `cat always_hold.txt | wc -l` atoms held; $silent_reswaps silent re-swaps of held atoms; $repulls re-pulls (not counted)"
endif
if( $test == 0 && ! $stuck && -e off_atoms.txt ) then
  if( ! -s off_atoms.txt ) then
    echo "arrived: maxD $arrive_maxd A"
    echo "all done."
    exit
  endif
endif
if( -s off_atoms.txt ) then
  echo "stuck: `cat off_atoms.txt | wc -l` atoms more than $arrive_tol A from the target:"
  head -n 20 off_atoms.txt
endif
if( $test == 0 ) then
  echo "not done: no relabel left, but not at the target."
else
  echo "not done. not sure what to do next."
endif

exit

########################################################################
#
#  notes and drafts below
#




foreach top ( `awk '{print ++n}' needed_swaps.txt` )

if( ! -s swaptop${top}_001.pdb ) continue
if( -e molprobify_swaptop${top}.log ) continue
echo $top

molprobify_runme.com swaptop${top}_001.pdb >&! molprobify_swaptop${top}.log &
sleep 1

end
wait








foreach top ( 0 `awk '{print ++n}' needed_swaps.txt` )

head -n $top needed_swaps.txt >! trial_swaps.txt

cat Hbonds.txt trial_swaps.txt $pdbfile |\
awk 'BEGIN{other["A"]="B";other["B"]="A"}\
  /^BOND/{hatoms[$2]=hatoms[$2]" "$3;next}\
  /^SWAP/{++swapconf[$2];\
    split($2,w,"_");atm=w[1];rid=w[2];split(hatoms[atm],a);\
    for(i in a){++swapconf[a[i]"_"rid];};next}\
  ! /^ATOM|^HETAT|^ANISO/{print;next}\
  {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
   gsub(" ","",atom);conf=substr($0,17,1);\
   pre=substr($0,1,16);post=substr($0,18);\
   atm=atom"_"chain resnum}\
  ! swapconf[atm]{print;next}\
  {print pre other[conf] post}' |\
cat >! swapped.pdb

if(-e swaptop${top}.log) continue

combine_pdbs_runme.com swapped.pdb $pdbfile outfile=swaptop${top}.pdb > /dev/null

$srun phenix.refine refme.mtz opts.eff swaptop${top}.pdb prefix=swaptop${top} >&! swaptop${top}.log &
sleep 1

end
wait

foreach top ( `awk '{print ++n}' needed_swaps.txt` )

if( ! -s swaptop${top}_001.pdb ) continue
if( -e molprobify_swaptop${top}.log ) continue
echo $top

molprobify_runme.com swaptop${top}_001.pdb >&! molprobify_swaptop${top}.log &
sleep 1

end
wait

sort -k2g *score.txt | head | awk '{print $1,$2,$3,$4,$(NF-2)}'

rm -f top_vs_maxd.txt
foreach top ( `awk '{print ++n}' needed_swaps.txt` )
  set maxd = `filter_pdb.awk -v skip=H,water good.pdb swaptop${top}_001.pdb | rmsd | awk '/^MAXD.all/{print }'`
  echo $top $maxd | tee -a top_vs_maxd.txt
end
set top = 23

head -n $top needed_swaps.txt >! sufficient_swaps.txt

awk '{print $2,sprintf("%.1f",$3)}' sufficient_swaps.txt |\
  awk -F "[_ ]" '! worst[$2]{worst[$2]=$3} {atoms[$2]=atoms[$2]" "$1}\
  END{for(r in atoms)print r,worst[r],atoms[r]}' |\
 sort -k2gr

foreach notswap ( `seq 1 $top` )

cat sufficient_swaps.txt | awk -v n=$notswap 'NR!=n{print}' >! trial_swaps.txt

cat Hbonds.txt trial_swaps.txt $pdbfile |\
awk 'BEGIN{other["A"]="B";other["B"]="A"}\
  /^BOND/{hatoms[$2]=hatoms[$2]" "$3;next}\
  /^SWAP/{++swapconf[$2];\
    split($2,w,"_");atm=w[1];rid=w[2];split(hatoms[atm],a);\
    for(i in a){++swapconf[a[i]"_"rid];};next}\
  ! /^ATOM|^HETAT|^ANISO/{print;next}\
  {atom=substr($0,12,5);chain=substr($0,22,1);resnum=substr($0,23,5)+0;\
   gsub(" ","",atom);conf=substr($0,17,1);\
   pre=substr($0,1,16);post=substr($0,18);\
   atm=atom"_"chain resnum}\
  ! swapconf[atm]{print;next}\
  {print pre other[conf] post}' |\
cat >! swapped.pdb

if(-e notswap${notswap}.log) continue

combine_pdbs_runme.com swapped.pdb $pdbfile outfile=notswap${notswap}.pdb

$srun phenix.refine refme.mtz opts.eff notswap${notswap}.pdb prefix=notswap${notswap} >&! notswap${notswap}.log &
sleep 1

end


rm -f notswap_vs_maxd.txt
foreach notswap ( `seq 1 $top` )
  set maxd = `filter_pdb.awk -v skip=H,water good_001.pdb notswap${notswap}_001.pdb | rmsd | awk '/^MAXD.all/{print }'`
  echo $notswap $maxd | tee -a notswap_vs_maxd.txt
end












exit


set mtzfile = best10_001.mtz

phenix.find_peaks_holes best10_001.mtz  best10_001.pdb

egrep "^ATOM" peaks_holes.pdb | head -n 1 >! atom.pdb

echo "labin F1=F0 PHI=PHI" | fft hklin refme.mtz mapout right.map

mtzdmp $mtzfile |\
  awk '/OVERALL/,/No. of refl/' |\
  awk 'NF<3{next} $(NF-1)=="F"{F=$NF} $(NF-1)=="P"{print F"|"$NF}' |\
cat >! labels.txt

foreach pair ( `cat labels.txt` )
  set F = `echo $pair | awk -F "|" '{print $1}'`
  set PHI = `echo $pair | awk -F "|" '{print $2}'`
  echo $F $PHI
  echo "labin F1=$F PHI=$PHI" | fft hklin $mtzfile mapout ${F}.map
end

 
cad hklin1 refme.mtz hklin2 $mtzfile hklout diff.mtz << EOF
labin file 1 E1=F0 E2=PHI
labin file 2 E1=F-model E2=PHIF-model
EOF

echo "labin F1=F0 PHI=PHI F2=F-model PHI2=PHIF-model" | fft hklin diff.mtz mapout diff.map


foreach prefix ( `ls -1rt *.mtz | awk -F "." '{print $1}'` )
if(! -e ${prefix}.pdb) continue
mmtbx.find_peaks_holes ${prefix}.pdb ${prefix}.mtz |& tee ${prefix}_peakhole.log &
end




