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

cat << EOF
pdbfile = $pdbfile
weight_snap = $weight_snap
fast = $fast
swapdown = $swapdown
gt_water = $gt_water
debug = $debug
EOF
which phenix.refine


set pdir = /home/jamesh/projects/build_challenge/1aho/
set mypath = `dirname $0`
setenv PATH ${mypath}:/home/jamesh/projects/amber.old/1aho_refine:${pdir}:${PATH}
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
    phenix.refine opts.eff refme.mtz refme.pdb prefix=optref main.number_of_macro_cycles=$cyc >! phenix1.log
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

# atoms whose swap required holdme.eff to secure (one-shot list);
# any of those that reappear in needed_swaps.txt get promoted to always_hold
# and are then swap-corrected + held in every subsequent refinement
if(! -e needs_holding.txt) then
   echo -n "" >! needs_holding.txt
endif
if(! -e always_hold.txt) then
   echo -n "" >! always_hold.txt
endif

set struggling = 0
set movetype = one_atom
set prev_movetype = $movetype
while ( -s needed_swaps.txt )

cp needed_swaps.txt needed_swaps_${swaps}.txt 

head -n 1 needed_swaps.txt | tee this_swap.txt
set atomid = `awk '{print $2;exit}' needed_swaps.txt`
set atom = `echo $atomid | awk -F "_" '{print $1}'`
set res = `echo $atomid | awk -F "_" '{print $2}'`
set resnum = `echo $res | awk '{print substr($1,2)}'`
set chain = `echo $res | awk '{print substr($1,1,1)}'`
set rigid = `echo $atom | awk '{print ( $1 ~ /[0-9]$/ || /CZ/ || /OH/ )}'`
set partner = `awk '/^SSBOND/{print $4 $5,$7 $8}' good.pdb |  awk '{print ;print $2,$1}' | egrep "^$res" | awk '{print $2}' | sort -u`

set this_atom = `awk '{print $2;exit}' needed_swaps.txt | sort`
set prev_atom = `awk '{print $2;exit}' needed_swaps_${prev}.txt | sort`
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
   set movetype = one_atom
   set prev_movetype = one_atom
endif

if( "$movetype" == "rigid_side" ) then
   echo "swapping rigid side chain of $res"
   cat needed_swaps.txt |\
   awk 'NR==1{split($2,w,"_");res=w[2]}\
     $2~/[1-9]_/ && $2~ "_"res"$"{print}' | tee this_swap.txt
endif
if( "$movetype" == "whole_res" ) then
    echo "swapping whole residue $res" 
    egrep "_${res} " needed_swaps.txt | tee this_swap.txt
endif
if( "$movetype" == "SSbond" ) then
    echo "swapping whole disulfide: $res - $partner" 
    egrep "_${res} |_${partner} " needed_swaps.txt | tee this_swap.txt
endif
if( "$movetype" == "top2" ) then
    echo "swapping top 2 entries" 
    head -n 2 needed_swaps.txt | tee this_swap.txt
endif
if( "$movetype" == "top3" ) then
    echo "swapping top 3 entries" 
    head -n 3 needed_swaps.txt | tee this_swap.txt
endif
if( "$movetype" == "all" ) then
    echo "swapping all entries"
    cat needed_swaps.txt | tee this_swap.txt
endif

# augment this_swap.txt with any always_hold atoms that are currently mis-labelled
# (i.e. they appear in needed_swaps.txt) so they get swap-corrected this iteration
if( -s always_hold.txt ) then
   awk '/^SWAP/{print $2}' this_swap.txt >! _curr_swap_set.txt
   awk 'NR==FNR{seen[$1]=1;next} ! seen[$1]' _curr_swap_set.txt always_hold.txt >! _ah_to_add.txt
   if( -s _ah_to_add.txt ) then
      awk 'NR==FNR{want[$1]=1;next} /^SWAP/ && want[$2]' _ah_to_add.txt needed_swaps.txt |\
         tee -a this_swap.txt | awk '{print "always_hold re-swap:",$0}'
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
compare_to_best_runme.com swapme.pdb prerefine=0 gt_water=$gt_water >&! compare_to_best_${swaps}.log
# updates lots of files

# any atom that previously needed hold and now reappears in needed_swaps.txt
# is promoted to always_hold (permanent swap-correction + reference restraint)
if( -s needs_holding.txt && -s needed_swaps.txt ) then
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

set test = `cat needed_swaps.txt | wc -l`
echo "$test swaps remain"
echo "$productive_steps productive moves in total"
if( $test == 0 ) then
  echo "all done."
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




