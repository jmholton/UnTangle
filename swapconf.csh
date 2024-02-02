#! /bin/tcsh -f
#
#   swap conformer letters of commanded residue range
#   optional: "side" or "main" or an atom name. I.E. "O" for just carbonyl oxygen
#
#
set pdbfile = refmacout.pdb
set outfile = swapped.pdb

set chain = "A"
set res1 = ""
set resn = ""
set opt = ""
# make chain letter same as conf letter
set chainconf = auto

set debug = 0
set tempfile = /dev/shm/${USER}/temp_$$_

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
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
      if("$key" == "debug") set debug = "$Val"
    else
      # no equal sign
      if("$arg" =~ *.pdb) then
          set pdbfile = "$Arg"
          continue
      endif

      # break up symbols
      set suf = "$arg"
      set ch  = ""
      if( "$Arg" =~ [A-Za-z][0-9]*) then
          set ch = `echo $Arg | awk '{print substr($0,1,1)}'`
          set suf = `echo $Arg | awk '{print substr($0,2)}'`
          set int = `echo $suf | awk '{print $1+0}'`
      endif

      if("$int" == "$suf") then
          if("$ch" != "") set chain = "$ch"
          if("$res1" == "") then
              set res1 = $int
              continue
          else
              set resn = $int
              continue
          endif
      endif
    endif
    # one-word-ers
    if("$arg" == "debug") set debug = 1
    if("$arg" == "main") set opt = "main"
    if("$arg" == "side") set opt = "side"
    if("$arg" == "rota") set opt = "rota"
    if("$arg" == "rotb") set opt = "rotb"
    if("$arg" == "all") set opt = "all"
    if("$Arg" =~ [CNOSH]*) set opt = "$Arg"
end

# figure out where to put temp files
if(! -w /dev/shm/ && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_sc_
if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_sc_
set tmpdir = `dirname $tempfile`
if(! -e "$tmpdir") mkdir -p "$tmpdir"

# shorthand
set t = ${tempfile}

# sanitize the range
if("$resn" == "") set resn = 0
if("$resn" == "0") set resn = $res1
if("$resn" == "") set resn = $res1

if(! -e "$pdbfile") then
   set BAD = "file $pdbfile does not exist\nusage: $0 model.pdb A1-15 CA\n to swap CA conformers in residues numbered 1 to 15"
   goto exit
endif

# work with a copy
egrep -v "^ANISO" $pdbfile >! ${t}old.pdb

# decide if we should swap chain as well as conf letter
if("$chainconf" == "auto") then
    awk '/^ATOM/{print substr($0,17,1),substr($0,22,1)}' ${t}old.pdb |\
    sort -u | sort >! ${t}confchains.txt
    set nchains = `awk '{print $2}' ${t}confchains.txt | sort -u | wc -l`
    set nconfs  = `awk '{print $1}' ${t}confchains.txt | sort -u | wc -l`
    set ncombo  = `cat ${t}confchains.txt | sort -u | wc -l`
    if( $ncombo == 2 && $nchains == 2 && $nconfs == 2 ) then
       set chainconf = 1
    else
       set chainconf = 0
    endif
endif

# now do the actual letter swaps
echo "swapping $opt atoms in chain $chain from $res1 to $resn"
echo "$res1 $resn $opt $chain" |\
cat - $pdbfile |\
awk 'NR==1{for(i=$1;i<=$2;++i)++sel[i];flag=$3;chain=$4;next}\
 /^ANISO/{next}\
 ! /^ATOM|^HETAT/{print;next}\
 {typ=substr($0,18,3);atm=substr($0,12,5);\
  pre=substr($0,1,16);post=substr($0,18);\
  c=substr($0,17,1);ch=substr($0,22,1);\
  resnum=substr($0,23,5)+0;\
  atom=atm;gsub(" ","",atm);\
  side=1}\
 {nc=" "} c=="A"{nc="B"} c=="B"{nc="A"}\
 flag=="all"{print pre nc post;next}\
 typ=="PRO" || atom~/  N  |  C  |  CA |  O  /{side=0}\
 ch!=chain{print;next}\
 ! sel[resnum]{print;next}\
 flag=="main" && side{print;next}\
 flag=="side" && ! side{print;next}\
 flag~/^[CNOSH]/ && flag != atm{print;next}\
 flag!~/^rot[ab]/{print pre nc post;next}\
 flag=="rota" && c!="A"{print;next}\
 flag=="rotb" && c!="B"{print;next}\
  typ=="HIS"{gsub(" CD2"," XXX");gsub(" ND1"," CD2");gsub(" XXX"," ND1");}\
  typ=="HIS"{gsub(" CE1"," XXX");gsub(" NE2"," CE1");gsub(" XXX"," NE2");}\
  typ=="ASP"{gsub(" OD1"," XXX");gsub(" OD2"," OD1");gsub(" XXX"," OD2");}\
  typ=="ASN"{gsub(" OD1"," XXX");gsub(" ND2"," OD1");gsub(" XXX"," ND2");}\
  typ=="GLU"{gsub(" OE1"," XXX");gsub(" OE2"," OE1");gsub(" XXX"," OE2");}\
  typ=="GLN"{gsub(" OE1"," XXX");gsub(" NE2"," OE1");gsub(" XXX"," NE2");}\
  {Ee=substr($0,14,1)}\
  typ~/HIS|ASN|GLN/ && Ee=="N"{gsub("   C","   N")}\
  typ~/HIS|ASN|GLN/ && Ee=="C"{gsub("   N","   C")}\
  typ=="TYR" || typ=="PHE"{gsub(" CD1"," XXX");gsub(" CD2"," CD1");gsub(" XXX"," CD2");}\
  typ=="TYR" || typ=="PHE"{gsub(" CE1"," XXX");gsub(" CE2"," CE1");gsub(" XXX"," CE2");}\
  typ=="TYR" || typ=="PHE"{gsub(" HD1"," XXX");gsub(" HD2"," HD1");gsub(" XXX"," HD2");}\
  typ=="TYR" || typ=="PHE"{gsub(" HE1"," XXX");gsub(" HE2"," HE1");gsub(" XXX"," HE2");}\
  {print}' |\
cat >! ${t}new.pdb

# and now re-organize the PDB file
if( $chainconf ) then
  foreach conf ( A B )
      echo $conf |\
      cat - ${t}new.pdb |\
      awk 'NR==1{c=$1;next}\
      /^CRYST/{print} ! /^ATOM|^HETAT/{next}\
      substr($0,17,1)==c{print substr($0,1,16),substr($0,18)}' |\
      cat >! ${t}conf${conf}.pdb
  end
  egrep "^CRYST" ${t}confA.pdb >! ${t}both.pdb
  foreach conf ( A B )
    echo $conf $chainconf |\
    cat - ${tempfile}conf${conf}.pdb |\
    awk 'NR==1{c=$1;chainconf=$2;next}\
       /^ATOM|^HETAT/{pre=substr($0,12,5);mid=substr($0,18,4);ch=substr($0,22,1);post=substr($0,23);\
       if(chainconf)ch=c;\
      printf("ATOM%7d%s%s%s%s%s\n",++n,pre,c,mid,ch,post)}' |\
    sort -k1.23,1.27g -k2,3g >> ${t}both.pdb
    rm -f ${t}conf${conf}.pdb
  end
  mv ${t}both.pdb ${t}new.pdb
endif
mv ${t}new.pdb $outfile
diff  ${t}old.pdb $outfile
ls -l $outfile

exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

if( ! $debug && "$tempfile" != "" && "$tempfile" != "./") then
   rm -f ${t}*
endif


exit







