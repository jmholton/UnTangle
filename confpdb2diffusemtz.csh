#! /bin/tcsh -f
#
#  calculate diffuse scatter due to two alt conformers, assuming every ASU is independent
#
#
set pdbfile = ""
set reso = 1.0
set CELL = ()
set SG = ""
set cellmult = 4

set ksol = 0.34
set Bsol = 40

set outfile = sqrtIdiffuse.mtz

set tempfile = /dev/shm/${USER}/temp_cp2d_$$_
#set tempfile = ./tempfile_Bud_
mkdir -p /dev/shm/${USER}
mkdir -p ${CCP4_SCR}
if( ! -e /dev/shm/${USER}) set tempfile = ./tempfile_cp2d_$$_

set logfile = details.log

echo "command-line arguments: $* "

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
      if("$key" == "output" || "$key" == "outmtz") set outfile = "$Val"
      if("$key" == "k_sol" ) set ksol = "$Val"
      if("$key" == "b_sol" ) set Bsol = "$Val"
      if("$key" == "mult" ) set cellmult = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) set pdbfile = "$Arg"
      if("$Arg" =~ *.mtz ) set mtzfile = "$Arg"
    endif
    if("$arg" == "debug") set debug = "1"
end

set t = "$tempfile"

if(! -e "$pdbfile") then
    set BAD = "pdbfile $pdbfile does not exist."
    goto exit
endif

cat << EOF
pdbfile = $pdbfile
outfile = $outfile
cellmult = $cellmult
tempfile = $tempfile
EOF

set pdbSG = `awk '/^CRYST/{print substr($0,56,12)}' $pdbfile | head -1`
set SG = `awk -v pdbSG="$pdbSG" -F "[\047]" 'pdbSG==$2{print;exit}' ${CLIBD}/symop.lib | awk '{print $4}'`
if("$SG" == "") then
    set SG = `echo $pdbSG | awk '{gsub(" ","");print}'`
    set SG = `awk -v SG=$SG '$4 == SG && $1 < 500 {print $4}' $CLIBD/symop.lib | head -1`
endif
if("$SG" == "") then
    set SG = `echo $pdbSG | awk '{gsub(" ","");print}'`
endif
# may need to be more clever here
set SG = `echo $SG | awk '{gsub("R","H"); print}'`
set pdbCELL = `awk '/^CRYST1/{print $2,$3,$4,$5,$6,$7}' $pdbfile`

if( $?CELL != 6 ) set CELL = ( $pdbCELL )

echo "getting symops for $SG"
awk -v SG=$SG '$4==SG{getline;while( /X/ ){print;getline}}' ${CLIBD}/symop.lib >! ${t}symops.txt


set confs = `awk '/^ATOM|^HETAT/{print substr($0,17,1)}' $pdbfile | sort -u | sort `
echo "found conformers: $confs"
set ops = `awk '{print NR}' ${t}symops.txt`


set bigcell = `echo $CELL $cellmult | awk '{print $NF*$1,$NF*$2,$NF*$3,$4,$5,$6}'`
echo "big cell: $bigcell"


foreach conf ( $confs )

echo $conf |\
cat - $pdbfile |\
awk 'NR==1{conf=$1;next}\
  ! /^ATOM|^HETAT/{next}\
  {c=substr($0,17,1)}\
  c==conf{print}' >! ${t}asu_${conf}.pdb

foreach op ( $ops )

set xyzop = `head -n $op ${t}symops.txt | tail -n 1`

echo "applying $xyzop in small cell to conf $conf"
pdbset xyzin ${t}asu_${conf}.pdb xyzout ${t}symgen_${conf}_${op}.pdb << EOF >> $logfile
CELL $CELL
symgen $xyzop
EOF

pdbset xyzin ${t}symgen_${conf}_${op}.pdb xyzout ${t}sfallme.pdb << EOF >> $logfile
cell $bigcell
SPACE 1
EOF

echo "fmodel"
rm -f ${t}symgen_${conf}_${op}.mtz
phenix.fmodel high_resolution=$reso ${t}sfallme.pdb \
   k_sol=$ksol b_sol=$Bsol \
   output.file_name=${t}symgen_${conf}_${op}.mtz >> $logfile

end
end

foreach op ( $ops )

echo "calculating FA^2+FB^2 vs (FA+FB)^2 for op $op"

rm -f ${t}diffuse_${op}.mtz
sftools << EOF >> $logfile
read ${t}symgen_${confs[1]}_${op}.mtz
read ${t}symgen_${confs[2]}_${op}.mtz
set labels
F1
P1
F2
P2
calc ( COL Fsum PHIsum ) = ( COL F1 P1 ) ( COL F2 P2 ) +
calc COL Favg = COL Fsum 2 /
calc COL Favgsq = COL Favg COL Favg *
calc COL I1 = COL F1 COL F1 *
calc COL I2 = COL F2 COL F2 *
calc COL Isum = COL I1 COL I2 +
calc COL Iavg = COL Isum 2 /
calc COL Idiff = COL Iavg COL Favgsq -
write ${t}diffuse_${op}.mtz col Idiff Favg Iavg
quit
y
EOF

end

echo "adding up..."
rm -f ${t}diffuse.mtz
foreach op ( $ops )

  echo $op
  if(! -e ${t}diffuse.mtz) then
     cp ${t}diffuse_${op}.mtz ${t}diffuse.mtz
     continue
  endif

rm -f ${t}new.mtz
sftools << EOF >> $logfile
read ${t}diffuse.mtz col Idiff
read ${t}diffuse_${op}.mtz col Idiff
set labels
I1
I2
calc COL Idiff = COL I1 COL I2 +
write ${t}new.mtz col Idiff
quit
y
EOF
mv ${t}new.mtz ${t}diffuse.mtz

end

echo "taking sqrt(Idiff)"
rm -f ${t}sqrt.mtz
sftools << EOF >> $logfile
read ${t}diffuse.mtz col Idiff
calc F COL sqrtIdiff = COL Idiff 0.5 **
write ${t}sqrt.mtz col sqrtIdiff
quit
y
EOF

mv ${t}sqrt.mtz ${outfile}

exit:

if( "${tempfile}" != "" && "${tempfile}" != "./" ) then
  rm -f ${tempfile}* >& /dev/null
endif

if($?BAD) then
    echo "ERROR $BAD"
    exit 9
endif


ls -l $outfile

exit



