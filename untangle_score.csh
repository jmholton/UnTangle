#! /bin/tcsh -f
#
#  throw every validation we can at it                          -James Holton 1-18-24
#
#  assemble an overall geometry score
#  "energy" = (deviate/sigma)**2
#
#  fraction of Gaussian-random events observed with at least one event above x=deviate/sigma 
#  in any of Nrep trials is: 
#  fracabove = 1-erf(deviate/sigma/sqrt(2))**Nreps
#
#   i.e. you expect to see > 2 sigma deviates ~5% of the time,
#   unless you look at 10 at a time, in which case 37% of the time you see at least one > 2 sigma deviate
#
#   so, probability that a x-sigma deviate is not noise, given N samples is:
#   Pnotnoise(deviate,N) = erf(abs(deviate)/sigma/sqrt(2))**N
#
#   final score will be sum of average square sigma deviations for each validation metric,
#   plus the sum of:
#            Pnotnoise*(worstdeviate/sigma)^2 
#    over worst outlier of each validation metric
#    to avoid overwhelming all other considerations, these outlier quantities are softened
#    so that energy values above 10 are substitued with 10+log(energy)-log(10)
#    also, the worst-clash energy is multiplied by the number of clashes
#
#   sigma level where probability is 0.5 is given by:
#      sigma_P50 = sqrt(2)*inverf((0.5)**(1./exponent))
#      which is reasonably approximated by:
#      a2 = -0.0341977; a1 = 0.852225 ; a0 = 1.10939
#      sigma_P50(log10e) = a2*log10e**2+a1*log10e+a0
#      log10e(exponent) = log(exponent)/log(10)
#      softer scoring function for large exponents is:
#      Pnotnoise = 1-2**((deviate/sigma_P50)**5) , naturally clipped at [-1:1]
#
#   for rama, rota etc.
#   convert probability/frequency back to sigma with: deviate/sigma = inverf(1-probOK)*sqrt(2)
#
#   cbetadev - use 0.05 A as "sigma"
#
#   nonbonds:  use Leonard-Jones to convert to energy [-1:inf], but dont let "worst" or avg be negative
#
#   omega twist: energy=((sin(omega)/0.07)^2+(1+cos(omega))^10)/(proxPRO*2+1)
#      where proxPRO means a neighboring residue is proline
#
#   allow "override" file    
#

set pdbfile = ""
set mtzfile = ""
set outprefix = "-"

set rstfile = ""
set topfile = ""

set tempfile = /dev/shm/${USER}/tempfile_$$_

set overridefile = ""

# for printing
set modulo = 9999
# ignore hydrogen bumps in ARG residues: ARGH
# ignore hydrogen bond lengths: HBON
# ignore water hydrogens: HOHH
set ignore = ARGH
# sigma deviate to include in refmac/phenix fudge files
set sigma_fudge = 3
# flag to write out restraint files for reducing outliers
set writefudge = 0

# flag to debug things
set debug = 0

foreach arg ( $* )
    set Key = `echo $arg | awk -F "=" '{print $1}'`
    set Val = `echo $arg | awk -F "=" '{print $2}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set val = `echo $Val | awk '{print tolower($1)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
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
      if("$key" == "outprefix") set outprefix = "$Val"
      if("$key" == "topfile") set topfile = "$Val"
      if("$key" =~ override*) set overridefile = "$Val"
      if("$key" == "ignarg") set ignore = ( $ignore ARGH )
      if("$key" == "keepgeo") set keepgeo = "$Val"
    else
      # no equal sign
      if("$Key" =~ *.pdb) set pdbfile = "$Key"
      if("$Key" =~ *.rst7) set rstfile = "$Key"
      if("$key" == "pdbfile") set pdbfile = "$Val"
      if("$key" == "outprefix") set outprefix = "$Val"
      if("$key" == "topfile") set topfile = "$Val"
      if("$key" =~ override*) set overridefile = "$Val"
      if("$key" == "modulo") set modulo = "$Val"
      if("$key" == "ignarg") set ignore = ( $ignore ARGH )
      if("$key" == "keepgeo") set keepgeo = "$Val"
      if("$key" == "tempfile") set tempfile = "$Val"
      if("$key" == "debug") set debug = "$int"
    endif
    if("$arg" == "debug") set debug = "1"
end

# check for dependencies
set test = `gnuplot --version | grep "gnuplot " | grep -v "Command not found" |wc -l`
if( ! $test ) then
   set BAD = "need gnuplot installed"
   goto exit
endif
# check for dependencies
set test = `phenix.version | grep "PHENIX: " | grep -v "Command not found" |wc -l`
if( ! $test ) then
   set BAD = "need phenix installed"
   goto exit
endif

if( $debug && "$tempfile" =~ /dev/shm/* ) set tempfile = ./tempfile
if( "$tempfile" == "" ) set tempfile = ./tmp_$$_
set tmpdir = ${tempfile}_dir
mkdir -p $tmpdir
if(! -w "$tmpdir") then
  set tempfile = ./tmp_$$_
  set tmpdir = ${tempfile}_dir
  mkdir -p $tmpdir
endif
set t = $tempfile

if("$mtzfile" == "" && -e "$pdbfile") then
    set base = `basename $pdbfile .pdb`
    set mtzfile = ${base}.mtz
endif
if("$outprefix" == "-" && -e "$pdbfile") then
    set base = `basename $pdbfile .pdb`
    set outprefix = ${base}
endif
if("$outprefix" == "-" && -e "$rstfile") then
   set base = `basename $rstfile .rst7`
    set outprefix = ${base}
endif
if("$outprefix" == "-") then
    set BAD = "no files"
    goto exit
endif

if(-e "$overridefile") then
    echo "using overrides in $overridefile"
endif

# default to unknown shifts
echo "x moved n/d A n/d occ n/d B n/d n/d n/d d" >! ${t}movement.txt
# "moved" dxyz "A" docc "occ" dB "B" atomX atomO atomB "d"
set refmacdir = `dirname $pdbfile`
foreach pref ( ${refmacdir}/refmac ${refmacdir}/${outprefix} ./refmac ./$outprefix  )
  set refmac_shifts = ${pref}_shifts.txt
  #ls "$refmac_shifts"
  if(-e "$refmac_shifts" ) break
end
if(-e "$refmac_shifts") then
  echo "reading $refmac_shifts"
  tail -n 3 $refmac_shifts |\
  awk 'NF>5{print "x",$0,"d"}' |\
  tail -n 1 >! ${t}movement.txt
  # "x moved" dxyz "A" docc "occ" dB "B" atomX atomO atomB "d"
endif

echo -n "" >! ${t}Rstats.txt
foreach pref ( ${refmacdir}/refmac ${refmacdir}/${outprefix} ./refmac ./$outprefix  )
  set refmac_shifts = ${pref}_shifts.txt
  set refmac_Rplot = ${pref}_Rplot.txt
  #ls "$refmac_Rplot"
  if(-e "$refmac_Rplot" ) break
end
if(-e "$refmac_Rplot") then
  echo "reading $refmac_Rplot"
  tail -n 3 $refmac_Rplot |\
  awk 'NF>11{n=$1+0;getline;$13+=0;print $0,n,"R"}' |\
  tail -n 1 >! ${t}Rstats.txt
  # trial Rw Rf FOM LL LLf rmsBond Zbond rmsAngle zAngle rmsChiral function vdw   trial    "R"
  #   1   2  3  4   5  6   7       8     9        10     11         12       13   14       15
endif

set phenix_log = `ls -1rt phenix*.log |& tail -n 1`
set test = `awk '$2+0>0' ${t}Rstats.txt | wc -l`
if(! $test && -e "$phenix_log") then
  echo "reading $phenix_log"
  tac $phenix_log |\
  awk '/Final R-work/{Rw=$4*100;Rf=$7*100}\
    $1=="end:"{Rw=$2*100;Rf=$3*100;bn=$4;ang=$5;}\
    /figures of merit:/{FOM=$NF+0}\
    /MACRO_CYCLE/{print 0,Rw+0,Rf+0,FOM+0,"n/d n/d",bn,"n/d",ang,"n/d n/d n/d n/d",0,"R";exit}' |\
  tail -n 1 >! ${t}Rstats.txt
  # trial Rw Rf FOM LL LLf rmsBond Zbond rmsAngle zAngle rmsChiral function vdw   trial    "R"
  #   1   2  3  4   5  6   7       8     9        10     11         12       13   14       15
endif

if(! -s  ${t}Rstats.txt) then
  echo "0 n/d n/d n/d n/d n/d n/d n/d n/d n/d n/d n/d 0 R" >! ${t}Rstats.txt
  # trial Rw Rf FOM LL LLf rmsBond Zbond rmsAngle zAngle rmsChiral function vdw   trial    "R"
  #   1   2  3  4   5  6   7       8     9        10     11         12       13   14       15
endif

set R_Rfree = `awk '/  R VALUE       |  FREE R VALUE   /{print $NF*100}' $pdbfile`
if( $#R_Rfree == 2 ) then
    echo "using R/Rfree from $pdbfile"
    echo $R_Rfree |\
    cat - ${t}Rstats.txt |\
    awk 'NR==1{R=$1;Rf=$2;next} {$2=R;$3=Rf;print}' |\
    cat >! ${t}.txt
    mv ${t}.txt ${t}Rstats.txt
endif

if(-e "$rstfile" && ! -e "$pdbfile") then
    set pdbfile = ${t}rst.pdb
    if(! $?AMBERHOME) source /programs/amber/amber.csh
    foreach dir ( ./ ../ ../../ )
    foreach name ( xtal.prmtop start.top )
      if(! -e "$topfile") set topfile = ${dir}$name
    end
    end
    if(! -e "$topfile") set topfile = `ls -1rt *top | tail -n 1`
    if(! -e "$topfile") then
        set BAD = "cannot find topology file"
        goto exit
    endif
    echo "using parmtop: $topfile"
    cpptraj << EOF > ${t}cpptraj.log
parm $topfile
trajin $rstfile
outtraj ${t}.pdb pdb pdbv3 sg "P 1"
go
EOF
    amb2pdb.awk ${t}.pdb >! $pdbfile
endif
set test = `cat $pdbfile | wc -l`
if( $test == 0 ) then
    set BAD = "empty PDB file"
    goto exit
endif
set test = `awk '$NF=="XP"' $pdbfile | wc -l`
if( $test ) then
    echo "discarding amber extra-point atoms"
    awk '$NF!="XP"' $pdbfile >! ${t}noxp.pdb
    set pdbfile = ${t}noxp.pdb
endif


#echo "phipsichi"
#phipsichi.com $pdbfile >! ${outprefix}_phipsichi.log

echo "ramalyze"
( phenix.ramalyze $pdbfile >! ${outprefix}_ramalyze.log ) >& /dev/null
echo "rotalyze"
( phenix.rotalyze $pdbfile >! ${outprefix}_rotalyze.log ) >& /dev/null
echo "omegalyze"
( phenix.omegalyze $pdbfile >! ${outprefix}_omegalyze.log ) >& /dev/null
echo "cbetadev"
( phenix.cbetadev $pdbfile  >! ${outprefix}_cbetadev.log ) >& /dev/null
#echo "clashscore" - will get from coot output from molprobity
#( molprobity.clashscore $pdbfile keep_hydrogens=True verbose=True --overwrite >! ${outprefix}_clashscore.log ) >& /dev/null
echo "molprobity"
set pwd = `pwd`
cp $pdbfile $tmpdir/this.pdb
( cd $tmpdir ;\
phenix.molprobity flip_symmetric_amino_acids=True \
    outliers_only=False output.probe_dots=False \
    output.coot=True this.pdb ) >&! ${outprefix}_molprobity.log
cp ${tmpdir}/molprobity_coot.py ${t}molprobity_coot.py

echo "geometry"
phenix.geometry_minimization $pdbfile macro_cycles=0 \
  output_file_name_prefix=${t} >! ${outprefix}_geom.log
# logfile is "greatest hits" only

# convert into more parsable format
cat ${t}.geo |\
awk '/nonbonded pdb=/{key="NONBOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     obs=$1;ideal=$2;sigma=1;energy=lj(obs,ideal);}\
   /bond pdb=/{key="BOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$4;energy=$6}\
   /angle pdb=/{key="ANGLE";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$4;energy=$6}\
   /dihedral pdb=/{key="TORSION";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;split($0,w,"\"");id4=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$5;energy=$7}\
   /chirality pdb=/{key="CHIR";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;split($0,w,"\"");id4=w[2];\
     getline;getline;\
     ideal=$2;obs=$3;sigma=$5;energy=$7}\
   /plane pdb=/{key="PLANE";ideal=0;obs=$(NF-4);sigma=$(NF-3);\
     while(NF && ! /sigma/){split($0,w,"\"");id=w[2];\
       a=substr(id,1,4);f=substr(id,5,1);t=substr(id,6,4);c=substr(id,10,1);r=substr(id,11,5);\
       if(f==" ")f="_";if(c==" ")c="_";\
       nid=sprintf("%4s %s %4s %s %5s",a,f,t,c,r);\
       if(! /plane/){obs=$(NF-2);sigma=$(NF-1)};\
       energy=(obs/sigma)^2;\
       print key,energy+0,ideal-obs,obs,ideal,sigma+0,"|",nid;\
       getline;}\
     }\
   id1!=""{\
     a1=substr(id1,1,4);a2=substr(id2,1,4);a3=substr(id3,1,4);a4=substr(id4,1,4);\
     f1=substr(id1,5,1);f2=substr(id2,5,1);f3=substr(id3,5,1);f4=substr(id4,5,1);\
     t1=substr(id1,6,4);t2=substr(id2,6,4);t3=substr(id3,6,4);t4=substr(id4,6,4);\
     c1=substr(id1,10,1);c2=substr(id2,10,1);c3=substr(id3,10,1);c4=substr(id4,10,1);\
     r1=substr(id1,11,5);r2=substr(id2,11,5);r3=substr(id3,11,5);r4=substr(id4,11,5);\
     gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
     gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
     nid1=sprintf("%4s %s %4s %s %5s",a1,f1,t1,c1,r1);\
     nid2=sprintf("%4s %s %4s %s %5s",a2,f2,t2,c2,r2);\
     nid3=sprintf("%4s %s %4s %s %5s",a3,f3,t3,c3,r3);\
     nid4=sprintf("%4s %s %4s %s %5s",a4,f4,t4,c4,r4);\
     ids=nid1; if(nid2~/[0-9]/)ids=ids" - "nid2; if(nid3~/[0-9]/)ids=ids" - "nid3;if(nid4~/[0-9]/)ids=ids" - "nid4;\
     print key,energy+0,ideal-obs,obs,ideal,sigma+0,"|",ids;\
     id1=id2=id3=id4=""} \
   function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
 function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
sort -k2gr |\
cat >! ${t}_fullgeo.txt

echo $ignore |\
cat - ${t}_fullgeo.txt |\
awk 'NR==1{$0==toupper($0);\
  if(/ARG/)++ARGH;\
  if(/HOH/)++HOHH;\
  if(/HBON/)++HBON;\
  next}\
  ARGH && /^NONBON/ && $8~/^H/ && $14~/^H/ && $10=="ARG" && $16=="ARG"{next}\
  HOHH && ( $8~/^H/ || $14~/^H/ ) && ( $10=="HOH" || $16=="HOH" ){next}\
  HBON && /^BOND/ && ( $8~/^H/ || $14~/^H/ ){next}\
  {print}' |\
cat >! ${t}_filtgeo.txt

cat ${t}_filtgeo.txt |\
awk '{atoms="";for(i=8;$i!="";i+=6){a=$i"_"$(i+4);atoms=atoms" "a};\
  print $1,$2,$3,$4,$5,$6,$7,atoms}' |\
tee ${t}_geoH.txt |\
awk '! /^NONBON/{for(i=7;i<=NF;++i)if($i~/^H/)next} {print}' |\
cat >! ${t}_geo.txt
# key energy delta obs ideal sigma | atomids

# fullgeo.txt:
# key energy delta obs ideal sigma | atom conf typ chain resnum - atom conf typ chain resnum - ...

# plane residual is sum of (delta/sigma)^2 for all involved atoms
# bond/angle/tors/chir residual is (delta/sigma)^2
# plane atoms should therefore be separate!  as they are above.
#

# get sequence
awk -F "|" '{print $2}' ${t}_fullgeo.txt |\
 awk -F "-" '{for(i=1;i<=NF;++i)print $i}' |\
 awk '{print $3,$4,$5}' |\
sort -u |\
sort -k2,2 -k3g >! ${t}sequence.txt


# convert other logs to parsable forms
cat ${outprefix}_omegalyze.log |\
awk -v modulo=$modulo -F ":" 'BEGIN{RTD=45/atan2(1,1)}\
   ! /^SUMMARY|^resid/{om=$3/RTD;\
   n=(substr($0,3,4)-1)%modulo+1;\
   energy=(sin(om)/0.07)^2+(1+cos(om))^10;\
   print "OMEGA",energy,n,$0}' |\
sort -k2gr >! ${t}_omegalyze.txt
# OMEGA energy resnum%64 otherstuff

# convert all omegas separately, noting prolines, treat "sigma" as 4 deg
awk '/PRO/{print "isPRO",$3}' ${t}sequence.txt |\
cat - ${t}_fullgeo.txt |\
awk -v modulo=$modulo 'BEGIN{RTD=45/atan2(1,1)} \
  /isPRO/{isPRO[$2]=1;next}\
  /^TORS/ && $8=="CA" && $26=="CA"{n=($12-1)%modulo+1;om=$4/RTD;\
    proxPRO=isPRO[n-1]+isPRO[n+1];\
    energy=((sin(om)/0.07)^2+(1+cos(om))^10)/(proxPRO*2+1);\
    print "OMEGA",energy,n,proxPRO,"omega=",om*RTD}' |\
sort -k2gr >! ${t}_allomegas.txt


# this is in Angstrom, assume "sigma" of 0.05 A
cat ${outprefix}_cbetadev.log |\
 awk -v modulo=$modulo -F ":" '! /^SUMM|^pdb/{energy=($6/0.05)^2;\
   n=($5-1)%modulo+1;\
   print "CBETADEV",energy,n,$0}' |\
sort -k2gr >! ${t}_cbetadev.txt
# CBETADEV energy resnum%modulo otherstuff


# convert rama, rota logs to probabilities too
# these are "% probability the rotamer is OK" 
cat ${outprefix}_rotalyze.log |\
 awk -v modulo=$modulo -F ":" '! /^SUMM|^resid/ && NF>8 && $2+0>0{\
   n=(substr($0,3,4)-1)%modulo+1;\
   P=1-$3/100;\
   pP=sprintf("%.35g",P-1e-16);\
   arg="inverf("pP")**2";\
   if(sprintf("%g",P)+0==1)arg=99;\
   print "print "arg", \"",n,$0,P+0,"\""}' |\
gnuplot |&\
awk '{$1+=0;print "ROTA",$0}' |\
sort -k2gr >! ${t}_rotalyze.txt
# ROTA energy resnum%64  otherstuff  1-probOK
cat ${outprefix}_ramalyze.log  |\
 awk -F ":" '! /^SUMM|^resid/ && /^ A/{\
   n=(substr($0,3,4)-1)%64+1;\
   P=1-$2/100;\
   pP=sprintf("%.35g",P-1e-16);\
   arg="inverf("pP")**2";\
   if(sprintf("%g",P)+0==1)arg=99;\
   print "print "arg", \"",n,$0,P,"\""}' |\
gnuplot |&\
awk '{$1+=0;print "RAMA",$0}' |\
sort -k2gr >! ${t}_ramalyze.txt
# RAMA energy resnum%64  otherstuff  1-probOK


# also get bad clashes
tail -n 2 ${t}molprobity_coot.py |\
 awk 'BEGIN{RS="[\047]"} {print}' |\
 awk 'length()==16{atoms=atoms"|"$0}\
 /^,/ && NF>1{r0=3;gsub(",","");r=r0+$1;\
  print "CLASH",lj(r,r0),-$1,atoms;atoms=""} \
  function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
  function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
sort -k2gr >! ${outprefix}_clashes.txt
# CLASH ljenergy deltadist "clash |" atoms1 "|" atoms2
#rm -f molprobity_coot.py


# make potential override list
if( 1 ) then
   cat ${t}_fullgeo.txt |\
   awk '{key=$1;obs=$4;ideal=$5;sigma=$6;\
         split($0,w,"|");n=split(w[2],w);atms="";\
         for(i=1;i<n;i+=6)atms=atms" "w[i]"_"w[i+4];\
         gsub("^ ","",atms);\
         target=obs}\
        /^NONBOND/ && obs>ideal {next}\
     {print key,"OVERRIDE",target,sigma,"|",atms}' |\
   cat >! potential_overrides.txt
   # TYPE "OVERRIDE" newtarget sigma "|" atoms
endif

# apply any override to preferred values
if(-e "$overridefile" ) then
    echo "applying overrides from $overridefile"
    cp ${t}_geo.txt ${t}_geo_orig.txt
    cat $overridefile ${t}_geo_orig.txt |\
    awk '{key=$1;energy=$2;delta=$3;obs=$4;ideal=$5;sigma=$6;\
     split($0,w,"|");atms=substr(w[2],match(w[2],/[^ ]/));\
     ent=$1" "atms}\
     $2=="OVERRIDE"{i=++orides[ent];oride_ideal[ent,i]=$3;oride_sigma[ent,i]=$4;next}\
     orides[ent]{n=orides[ent];\
        for(i=1;i<=n;++i){\
          t=oride_ideal[ent,i];s=oride_sigma[ent,i];\
          if(t=="")t=ideal;\
          if(s=="")s=sigma;\
          dev=(t-obs);newE=(dev/s)^2+energy/1e6;\
          if(key=="NONBOND"){\
            if(t>obs)t=ideal;\
            newE=lj(obs,t)};\
          if(newE<energy){ideal=t;energy=newE;delta=dev;ideal=t;sigma=s};\
        }\
        print key,energy,delta,obs,ideal,sigma,"|",atms;next}\
     {print}\
     function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
     function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
     sort -k2gr >! ${t}_geo.txt


    cp ${t}_fullgeo.txt ${t}_fullgeo_orig.txt
    cat $overridefile ${t}_fullgeo_orig.txt |\
    awk '{key=$1;energy=$2;delta=$3;obs=$4;ideal=$5;sigma=$6;\
     split($0,w,"|");atoms=w[2];gsub("^ ","",atoms);\
     n=split(atoms,w);atms="";\
     for(i=1;i<n;i+=6)atms=atms" "w[i]"_"w[i+4];\
     gsub("^ ","",atms);\
     ent=$1" "atms;\
     if(debug) print "DEBUG1",ent,n;\
     }\
     $2=="OVERRIDE"{atms=atoms;\
     ent=$1" "atms;\
     i=++orides[ent];oride_ideal[ent,i]=$3;oride_sigma[ent,i]=$4;\
     if(debug) print "DEBUG2",ent;\
     next}\
     orides[ent]{n=orides[ent];\
        for(i=1;i<=n;++i){\
          t=oride_ideal[ent,i];s=oride_sigma[ent,i];\
          if(t=="")t=ideal;\
          if(s=="")s=sigma;\
          dev=(t-obs);newE=(dev/s)^2+energy/1e6;\
          if(key=="NONBOND"){\
            if(t>obs)t=ideal;\
            newE=lj(obs,t)};\
          if(newE<energy){\
            if(debug) print "DEBUG3 overriding",ent;\
            ideal=t;energy=newE;delta=dev;ideal=t;sigma=s;\
          };\
        }\
        print key,energy,delta,obs,ideal,sigma,"|",atoms;next}\
     {print}\
     function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2^(-1./6)/r)^12-(r0*2^(-1./6)/r)^6)}\
     function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
     sort -k2gr >! ${t}_fullgeo.txt



     # edit other tables to ignore overriden values
     egrep "^TORS" $overridefile |\
     awk '{atoms=substr($0,index($0,"|")+1);\
       split(atoms,a);split(a[1],w,"_");n=w[2];\
       gsub("[ 0-9]","",atoms);\
       gsub("_"," ",atoms);print n,atoms}' |\
     awk '{type="rama"}\
          $2!~/^CA$|^N$|^C$/{type="rota"}\
          $3!~/^CA$|^N$|^C$/{type="rota"}\
          $4!~/^CA$|^N$|^C$/{type="rota"}\
          $5!~/^CA$|^N$|^C$/{type="rota"}\
          $2=="CA" && $NF=="CA"{type="omega"}\
       {print $0,type}' |\
     cat >! ${t}_orides.txt
     # resnum  a1 a2 a3 a4  rota/omega
     # shouldnt actually be any rama restraints

     cp ${t}_rotalyze.txt ${t}_rotalyze_orig.txt
     awk '$NF=="rota"{print $1}' ${t}_orides.txt |\
     sort -u | sort -g |\
     cat - ${t}_rotalyze_orig.txt |\
     awk 'NF==1{++ignore[$1];next}\
        ignore[$3]{next} {print}' |\
     cat >! ${t}_rotalyze.txt
     # ROTA energy resnum%64  otherstuff  1-probOK

     cp ${t}_ramalyze.txt ${t}_ramalyze_orig.txt
     awk '$NF=="rama"{print $1}' ${t}_orides.txt |\
     sort -u | sort -g |\
     cat - ${t}_ramalyze_orig.txt |\
     awk 'NF==1{++ignore[$1];next}\
        ignore[$3]{next} {print}' |\
     cat >! ${t}_ramalyze.txt
     # RAMA energy resnum%64  otherstuff  1-probOK

     cp ${t}_omegalyze.txt ${t}_omegalyze_orig.txt
     awk '$NF=="omega"{print $1}' ${t}_orides.txt |\
     sort -u | sort -g |\
     cat - ${t}_omegalyze_orig.txt |\
     awk 'NF==1{++ignore[$1];next}\
        ignore[$3]{next} {print}' |\
     cat >! ${t}_omegalyze.txt

     cp ${t}_allomegas.txt ${t}_allomegas_orig.txt
     awk '$NF=="omega"{print $1}' ${t}_orides.txt |\
     sort -u | sort -g |\
     cat - ${t}_allomegas_orig.txt |\
     awk 'NF==1{++ignore[$1];next}\
        ignore[$3]{next} {print}' |\
     cat >! ${t}_allomegas.txt
     # OMEGA energy resnum%64 proximal2PRO otherstuff omega

     # override clashes? cbetadevs? 
endif
# bond, angle, tors : residual= (delta/sigma)^2
# planes: sum up (dev/sigma)^2 for all atoms involved





# now that all exceptions/overrides have been applied ...

# set variables for worst of each type

# find worst non-bond
awk '/^NONBON/ && $2>-0.5 && $4<$5{print}' ${t}_filtgeo.txt |\
sort -k2gr >! ${outprefix}_badnonbond.txt
# NONBOND ljenergy delta obs ideal | atomids 

# make sure at least one line
set test = `wc -l ${outprefix}_badnonbond.txt | awk '{print $1}'`
if( ! $test ) then
    awk '/^NONBON/{print;exit}' ${t}_filtgeo.txt |\
    cat >! ${outprefix}_badnonbond.txt
endif

# and worst one that counts
awk '/^NONBON/ && $2>-0.1 && $4<$5{print}' ${t}_geo.txt |\
sort -k2gr >! ${outprefix}_scorednonbond.txt
# NONBOND ljenergy delta obs ideal | atomids 


# make sure at least one line
set test = `wc -l ${outprefix}_scorednonbond.txt | awk '{print $1}'`
if( ! $test ) then
    awk '/^NONBON/{print;exit}' ${t}_geo.txt |\
    cat >! ${outprefix}_scorednonbond.txt
endif

cat ${outprefix}_scorednonbond.txt >! ${t}_worstnonbond.txt
set worst_nonbond = `awk '{$1="";print;exit}'  ${t}_worstnonbond.txt`
# ljenergy delta obs ideal | atomids 


# energy delta obs ideal sigma | atoms
egrep "^BOND" ${t}_geo.txt | sort -k2gr >! ${t}_worstbond.txt
set worst_bond = `awk '{$1="";print;exit}' ${t}_worstbond.txt`

egrep "^ANGLE" ${t}_geo.txt | sort -k2gr >! ${t}_worstangle.txt
set worst_angle = `awk '{$1="";print;exit}'  ${t}_worstangle.txt`

egrep "^TORS" ${t}_geo.txt | sort -k2gr >! ${t}_worsttorsion.txt
set worst_torsion = `awk '{$1="";print;exit}'  ${t}_worsttorsion.txt`

egrep "^CHIR" ${t}_geo.txt | sort -k2gr >! ${t}_worstchiral.txt
set worst_chiral = `awk '{$1="";print;exit}'  ${t}_worstchiral.txt`

egrep "^PLANE" ${t}_geo.txt | sort -k2gr >! ${t}_worstplane.txt
set worst_plane = `awk '{$1="";print;exit}'  ${t}_worstplane.txt`

# worst omega outlier counts extra
set worst_omega = `awk '{$1="";print;exit}'  ${t}_allomegas.txt`


# worst outliers in rota rama and CBetadev
set worst_rota = `awk '{$1="";print;exit}' ${t}_rotalyze.txt`
set worst_rama = `awk '{$1="";print;exit}' ${t}_ramalyze.txt`
set worst_CB = `awk '{$1="";print;exit}' ${t}_cbetadev.txt`
# key energy resnum%64  otherstuff 
awk '/^CBET/{print $0,"| CB_"$3}' ${t}_cbetadev.txt >! ${t}_worstcbetadev.txt


# CLASH ljenergy deltadist "|" atoms1 "|" atoms2
# multiply worst clash by number of clashes
set nclashes = `cat ${outprefix}_clashes.txt | wc -l`
set worst_clash = `awk -v n=$nclashes '{$1=n"x";print sqrt(1+$2*n)^2;print;exit}' ${outprefix}_clashes.txt`
if("$worst_clash" == "") set worst_clash = "no"
# nclash*worst CLASH worstenergy nclash delta | atomid1 | atomid2


# now do R factor - assume "sigma" is 3% ?






# if missing
if("$worst_nonbond" == "") set worst_nonbond = "n/d"
if("$worst_clash" == "") set worst_clash = "no"
if("$worst_bond" == "") set worst_bond = "n/d"
if("$worst_angle" == "") set worst_angle = "n/d"
if("$worst_torsion" == "") set worst_torsion = "n/d"
if("$worst_chiral" == "") set worst_chiral = "n/d"
if("$worst_plane" == "") set worst_plane = "n/d"

if("$worst_omega" == "") set worst_omega = "n/d"
if("$worst_rota" == "") set worst_rota = "n/d"
if("$worst_rama" == "") set worst_rama = "n/d"
if("$worst_CB" == "") set worst_CB = "n/d"



# in addition to worsts, use average as part of score

# average "energy" residuals
set avg_bond    = `awk '/^BOND/ {++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_nonbond = `awk '/^NONBOND/ && $2>0{++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_full_nonbond = `awk '/^NONBOND/ {++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_angle   = `awk '/^ANGLE/{++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_torsion = `awk '/^TORS/ {++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_chiral  = `awk '/^CHIR/ {++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`
set avg_plane   = `awk '/^PLANE/{++n;sum+=$2} END{print sum/n,n}' ${t}_geo.txt`

set avg_omega = `awk '{sum+=$2;++n} END{if(n) print sum/n,n}' ${t}_allomegas.txt`
set avg_rama  = `awk '{sum+=$2;++n} END{if(n) print sum/n,n}' ${t}_ramalyze.txt`
set avg_rota  = `awk '{sum+=$2;++n} END{if(n) print sum/n,n}' ${t}_rotalyze.txt`
set avg_CB    = `awk '{sum+=$2;++n} END{if(n) print sum/n,n}' ${t}_cbetadev.txt`

set avg_clash = `awk '{sum+=$2;++n} END{if(n) print sum/n,n}' ${outprefix}_clashes.txt`

# fill in missing stats
if( "$avg_bond"    == "" ) set avg_bond = "n/d 0"
if( "$avg_nonbond" == "" ) set avg_nonbond = "n/d 0"
if( "$avg_full_nonbond" == "" ) set avg_full_nonbond = "n/d 0"
if( "$avg_angle"   == "" ) set avg_angle = "n/d 0"
if( "$avg_torsion" == "" ) set avg_torsion = "n/d 0"
if( "$avg_chiral"  == "" ) set avg_chiral = "n/d 0"
if( "$avg_plane"   == "" ) set avg_plane = "n/d 0"

if( "$avg_omega" == "" ) set avg_omega = "n/d 0"
if( "$avg_rota" == "" )  set avg_rota = "n/d 0"
if( "$avg_rama" == "" )  set avg_rama = "n/d 0"
if( "$avg_CB" == "" )    set avg_CB = "n/d 0"

if( "$avg_clash" == "" ) set avg_clash = "no 0"




cat << EOF |& tee ${t}_avgs.txt
$avg_bond    avg_bond    
$avg_nonbond avg_nonbond 
$avg_full_nonbond avg_full_nonbond 
$avg_clash   avg_clash   
$avg_angle   avg_angle   
$avg_torsion avg_torsion  
$avg_chiral  avg_chiral    
$avg_plane   avg_plane   
$avg_omega   avg_omega   
$avg_rota    avg_rota
$avg_rama    avg_rama    
$avg_CB      avg_CBdev
EOF

cat << EOF | tee ${t}_worsts.txt
$worst_bond     bond
$worst_nonbond  nonbond
$worst_clash    clash
$worst_angle    angle
$worst_torsion  torsion
$worst_chiral   chiral
$worst_plane    plane
$worst_omega    omega
$worst_rama     rama
$worst_rota     rota
$worst_CB       CBdev
EOF

#cat << EOF | tee ${t}_Rstats.txt
#EOF




# convert energy and frequency to probability its not noise (Pnn)
cat << EOF >! ${t}gnuplot.in
y0              = 1
a2              = -0.0192266
a1              = 0.751694
a0              = 1.12482
mx              = 0.21805
my              = 0.736621
# probabilty of not being noise
Pnn(deviate,N) = erf(abs(deviate)/sqrt(2))**N
# softer approximation
softPnna(deviate,N)= 1-2.0**clip((-abs(deviate/asigma_Pnn50(safelog(N)))**exponent(deviate,safelog(N))),1000)
# clip x to +/- y to range
clip(x,y) = (x>y?y:(x<=-y?-y:x))
# sigma deviates that have Pnn=0.5
sigma_Pnn50(N) = sqrt(2)*inverf((0.5)**(1./N))
asigma_Pnn50(logN) = (a2*logN**2+a1*logN+a0)
# exponent to use in approximation
exponent(deviate,logN) = my*logN+mx*deviate+y0
safelog(x) = (x<=0?0:log10(x))
s_thresh(N) = sqrt(2)*inverf((1e-30)**(1./N))
safePnn(deviate,N) = (deviate<s_thresh(N)?0:Pnn(deviate,N))

# chi square test
chisq(x,k) = (k<1?0:(k<3200?igamma(k/2.0,x/2.0):chisqhi(x,k)))
chisqhi(x,k) = erf((sqrt(x)-sqrt(k)))/2+0.5
EOF

cat ${t}_avgs.txt ${t}_worsts.txt |\
awk '{tag=$NF}\
   {$0=tolower($0);gsub("cbetadev","cbdev");gsub("torsion","tors");gsub("chiral","chir")}\
   {key=$NF;gsub("avg_","",key)}\
   {energy=$1} energy<0{energy=0} energy>400{energy=400}\
   tag~/^avg/{Nrep[key]=$2;\
   ssd=energy*Nrep[key];\
   print "print chisq("ssd","Nrep[key]"),\" Chisq ",tag,n"\"";next}\
   Nrep[key]==""{Nrep[key]=1}\
   {delta=sqrt(energy);\
    print "print softPnna("delta","Nrep[key]"),\" Pnn ",tag"\"";}' |\
cat >! ${t}print.gnuplot

cat ${t}gnuplot.in ${t}print.gnuplot |\
gnuplot >&! ${outprefix}_weights.txt

# some kind of sanity check here?

cat ${outprefix}_weights.txt ${t}_avgs.txt ${t}_worsts.txt |\
awk '$2~/[A-Z]/{weight[$NF]=$1;next}\
  {energy=$1+0}\
  energy<=0{energy=0}\
  energy>10{energy=10+log(energy)-log(10)}\
  /clash/{energy+=$2}\
  {print energy*weight[$NF],$0}' |\
sort -gr |\
cat >! ${t}_weighted.txt
# w_energy tag

# sum of all weighted deviates
set energy = `awk '{sum+=$1} END{print sum}' ${t}_weighted.txt`
echo "weighted energy (wE): $energy"


# save for later
awk '! /^NONBON/ && $2>4' ${t}_fullgeo.txt |\
cat >! ${outprefix}_worstgeo.txt

if( $?keepgeo ) then
  cat ${t}_fullgeo.txt |\
  awk '{atoms="";for(i=8;$i!="";i+=6){a=$i"_"$(i+4);atoms=atoms" "a};\
    print $1,$2,$3,$4,$5,$6,$7,atoms}' |\
  tee ${t}_geoH.txt |\
  awk '! /^NONBON/{for(i=7;i<=NF;++i)if($i~/^H/)next} {print}' |\
  cat >! ${outprefix}_geo.txt
  cp ${t}_fullgeo.txt ${outprefix}_fullgeo.txt
endif


if( $writefudge == 0 ) goto skipfudge
# now convert to refmac restraints

# translate hydrogens to bonded parents
cat ${t}_fullgeo.txt |\
awk '! /^BOND/{next}\
   $8~/^H/{print "PARENT",$8,$10,"-",$14}\
   $14~/^H/{print "PARENT",$14,$16,"-",$8}' |\
sort -u >! ${t}_parents.txt

# grab all non-bad atoms and hold in place
awk '$2>0{print $2+0.0001,"-"substr($0,index($0,"|")+1)}' ${t}_fullgeo.txt |\
awk -F "-" '{for(a=2;a<=NF;++a)print $1,"-"$a}' |\
awk '{while(gsub(/ $/,""));while(gsub("  "," "));print}' |\
cat ${t}_parents.txt - |\
awk '/^PARENT/{parent[$2,$3]=$5;next}\
  parent[$3,$5]{$3=parent[$3,$5]}\
  {print}' |\
awk -F "-" '$1>worst[$2]+0{worst[$2]=$1}\
   END{for(a in worst) print worst[a],a}' |\
sort -u |\
sort -g |\
awk '$1<0.1{$1=0.1}\
   $1<1{\
   a1=$2;f1=$3;t1=$4;c1=$5;r1=$6;\
   sigma=$1;\
   print "extern harmonic atinfo chain",c1,"resi",r1,"atom",a1,"alt",f1,"sigma",sigma;\
  }' |\
cat >! refmac_opts_holdme.txt



# make anti-clash restraints
egrep "^CLASH" ${outprefix}_clashes.txt |\
awk '{print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"|")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   a1=substr(id1,13,5);a2=substr(id2,13,5);\
   f1=substr(id1,8,1);f2=substr(id2,8,1);\
   t1=substr(id1,9,4);t2=substr(id2,9,4);\
   c1=substr(id1,2,1);c2=substr(id2,2,1);\
   r1=substr(id1,3,5);r2=substr(id2,3,5);\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   print "select",$3,a1,f1,t1,c1,r1,a2,f2,t2,c2,r2}' |\
cat ${t}_parents.txt - |\
awk '/^PARENT/{parent[$2,$3]=$5;next}\
  parent[$3,$5]{$3=parent[$3,$5]}\
  parent[$8,$10]{$8=parent[$8,$10]}\
  {print}' |\
egrep "^select" >! ${t}_selections.txt

cat ${t}_selections.txt ${t}_fullgeo.txt |\
awk '/^select/{del[$3,$4,$5,$6,$7]=del[$8,$9,$10,$11,$12]=$2;next}\
  ! /^NONBON/{next}\
  {nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
  }\
   del[a1,f1,t1,c1,r1] && del[a2,f2,t2,c2,r2]{\
   print ;\
   ideal=$5+del[a1,f1,t1,c1,r1];sigma=0.2;\
   print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999";\
   print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0"}' |\
awk '{gsub("altecode _","");gsub("chain _","");print}' |\
egrep "^extern|^#extern" >! refmac_opts_declash.txt


cat ${t}_selections.txt ${t}_fullgeo.txt |\
awk '/^select/{del[$3,$4,$5,$6,$7]=del[$8,$9,$10,$11,$12]=$2;next}\
  ! /^NONBON/{next}\
  {nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
  }\
   del[a1,f1,t1,c1,r1] && del[a2,f2,t2,c2,r2]{\
   print "#",$0;\
   ideal=$5+del[a1,f1,t1,c1,r1];sigma=0.2;\
   deadband=1;\
   print "    bond {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      distance_ideal =",ideal+deadband;\
   print "      sigma =",sigma;\
   print "      slack =",deadband/2;\
   print "    }"}' |\
awk '{gsub("and altid _","");gsub("and chain _","");print}' |\
cat >! phenix_declash.txt







# key ljenergy delta obs ideal 1 | atoms1 | atoms2
awk '/^NONBON/ && $2>1' ${t}_fullgeo.txt |\
cat ${t}_parents.txt - |\
awk '/^PARENT/{parent[$2,$3]=$5;next}\
  {print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   energy=$2;delta=$3;ideal=$5;\
   if(parent[a1,t1]){a1=parent[a1,t1];delta+=0.1}\
   if(parent[a2,t2]){a2=parent[a2,t2];delta+=0.1}\
   print "select",energy,delta,ideal,a1,f1,t1,c1,r1,a2,f2,t2,c2,r2}' |\
egrep "^select" >! ${t}_selections.txt
# parent atoms selected

# now look for distance between parent atoms in the non-bond list
# if its not there, make something up
cat ${t}_selections.txt ${t}_fullgeo.txt |\
awk '/^select/{\
     a1 =$5;f1 =$6;t1 =$7;c1 =$8;r1 =$9;\
     a2=$10;f2=$11;t2=$12;c2=$13;r2=$14;\
     key=a1" "f1" "c1" "r1" "a2" "f2" "c2" "r2;\
     oldE[key]=$2;\
     olddel[key]=$3;\
     oldideal[key]=$4;\
     ++selected[key];next}\
  ! /^NONBON/{next}\
  {nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   energy=$2;delta=$3;ideal=$5;\
   key=a1" "f1" "c1" "r1" "a2" "f2" "c2" "r2;\
   #print del[key]+0,selected[key],key;\
  }\
  selected[key]{\
   print ;\
   dd=del[key];\
   ideal+=dd;delta+=dd;\
   sigma=0.4/(delta+1.1);\
   if(sigma<0.001)sigma=0.001;\
   print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999";\
   print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0";\
   ++printed[key]}\
  END{for(key in selected){\
   if(selected[key] && ! printed[key]){\
     split(key,w);\
     a1=w[1];f1=w[2];c1=w[3];r1=w[4];\
     a2=w[5];f2=w[6];c2=w[7];r2=w[8];\
     sigma=0.4/(olddel[key]+1.1);\
     ideal=oldideal[key]+olddel[key]+0.7;\
     print selected[key],( ! selected[key] ),del[key],key;\
     print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
        "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
        "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999";\
     print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
        "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
        "value",ideal,"sigma",sigma,"type 0";\
   };\
  }}' |\
awk '{gsub("altecode _","");gsub("chain _","");print}' |\
egrep "^extern|^#extern" >! refmac_opts_debump.txt




# key energy delta obs ideal sigma | atoms
egrep "^BOND " ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1,i);\
   print "|"id1"|";\
   print "|"id2"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   ideal=$5;sigma=$6;\
   if(sigma<0.001)sigma=0.001;\
   if(a1~/^H/ || a2~/^H/){next}\
   print "extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0"}' |\
awk '{gsub("altecode _","");print}' |\
egrep "^extern" >! refmac_opts_fixbond.txt



# key energy delta obs ideal sigma | atoms
egrep "^ANGLE " ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   i=index($0,"|");\
   nxt=substr($0,i+1);\
   i=index(nxt,"-");\
   id1=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id2=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id3=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   print "|"id3"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   split(id3,w);a3=w[1];f3=w[2];t3=w[3];c3=w[4];r3=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
   gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
   ideal=$5;sigma=$6;\
   print "extern angle first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "next chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "next chain",c3,"resi",r3,"atom",a3,"altecode",f3,\
      "value",ideal,"sigma",sigma}' |\
awk '{gsub("altecode _","");print}' |\
egrep "^extern" >! refmac_opts_fixang.txt


# add cbetadev - tighten up angles around CA
# CBETADEV energy resnum%modulo otherstuff
egrep "^CBETA" ${t}_cbetadev.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   n=split($0,w,":");\
   f=w[2];c=w[4];r=w[5];\
   ideal=$5;sigma=$6;\
   print "extern angle first chain",c,"resi",r,"atom",a,"altecode",f,\
      "next chain",c,"resi",r,"atom",a,"altecode",f,\
      "next chain",c,"resi",r,"atom",a,"altecode",f,\
      "value",ideal,"sigma",sigma}' |\
awk '{gsub("altecode _","");print}' |\
egrep "^extern" >! refmac_opts_fixCBdev.txt
rm -f refmac_opts_fixCBdev.txt



# key energy delta obs ideal sigma | atoms
egrep "^TOR" ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   i=index($0,"|");\
   nxt=substr($0,i+1);\
   i=index(nxt,"-");\
   id1=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id2=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id3=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id4=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   print "|"id3"|";\
   print "|"id4"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   split(id3,w);a3=w[1];f3=w[2];t3=w[3];c3=w[4];r3=w[5];\
   split(id4,w);a4=w[1];f4=w[2];t4=w[3];c4=w[4];r4=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
   gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
   ideal=$5;sigma=$6;\
   print "extern torsion first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "next chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "next chain",c3,"resi",r3,"atom",a3,"altecode",f3,\
      "next chain",c4,"resi",r4,"atom",a4,"altecode",f4,\
      "value",ideal,"sigma",sigma,"period 0"}' |\
awk '{gsub("altecode _","");print}' |\
egrep "^extern" >! refmac_opts_fixtorsion.txt








# all the way back to geo file to do planes
cat ${t}.geo |\
awk '/plane pdb=/{sigma=$(NF-3)+0;rmsd=$(NF-1)+0;resid=$NF+0;n=0;\
     if(rmsd/sigma<sigcut)next;\
     printf("PLANE %s %s ",rmsd/sigma,sigma);\
     while(NF && ! /sigma/){split($0,w,"\"");id=w[2];\
       a=substr(id,1,4);f=substr(id,5,1);t=substr(id,6,4);c=substr(id,10,1);r=substr(id,11,5);\
       if(c==" ")c="_";if(f==" ")f="_";\
       atm=a;gsub("^ ","",atm);\
       if(atm !~ /^H/) printf("%s %s %s %s  ",a,f,c,r);\
       getline;\
     };\
     print "";}' |\
cat >! ${t}allplanes.txt

cat ${t}_worstplane.txt |\
awk -v s=$sigma_fudge '$2>s {print substr($0,index($0,"|")+1)}' |\
cat - ${t}allplanes.txt |\
awk 'NF==1 && $1 !~ /^H/{++sel[$1];next}\
  {thisplane=0;for(i=4;i<NF;i+=4){atm=$i"_"$(i+3);\
    if(sel[atm])++thisplane;}}\
  thisplane{print}' |\
cat >! ${t}selected_planes.txt

cat ${t}selected_planes.txt |\
awk '/^PLANE/ && NF>15{sigma=$3;\
   print;\
   printf("extern plane ");\
   for(i=4;i+4<NF;i+=4){\
     a=$i;f=$(i+1);c=$(i+2);r=$(i+3);\
     printf("next chain %s resi %s atom %s altecode %s ",c,r,a,f);}\
   print "sigma",sigma}' |\
awk '{gsub("altecode _","");gsub("chain _","");gsub("plane next","plane first");print}' |\
egrep "^extern" >! refmac_opts_fixplanes.txt


cat refmac_opts_fixplanes.txt |\
awk '/extern plane/{sigma=$NF;\
     print "    planarity {";\
     printf("      atom_selection = \"");\
     for(i=5;i<NF;i+=9){\
       c=$i;r=$(i+2);a=$(i+4);f=$(i+6);\
       printf("( name %s and resseq %s and chain %s and altid %s) ",a,r,c,f);\
       if(i<NF-9){printf("\n        or ")};\
     }\
     print "\"\n      sigma =",sigma;\
     print "    }";\
   }' |\
cat >! phenix_opts_fixplanes.txt


foreach suff ( declash debump )
 egrep "^#extern" refmac_opts_${suff}.txt |\
 awk -v deadband=1 '$2=="dist"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     v=$22;s=$24;\
   print "    #",c1,c2,a1,a2,r1,r2,f1,f2,v,s;\
   print "    bond {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      distance_ideal =",v+deadband;\
   print "      sigma =",s;\
   print "      slack =",deadband;\
   print "    }";}' |\
 cat >! phenix_opts_${suff}.txt
end

egrep "^extern" refmac_opts_fixang.txt |\
awk '$2=="angle"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     c3=$23;r3=$25;a3=$27;f3=$29;\
     v=$31;s=$33;\
   print "    #",c1,c2,c3,a1,a2,a3,r1,r2,r3,f1,f2,f3,v,s;\
   print "    angle {"\
   print "      action = *change";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      atom_selection_3 = \"name",a3,"and resseq",r3,"and chain",c3,"and altid",f3 "\"";\
   print "      angle_ideal =",v;\
   print "      sigma =",s;\
   print "    }";\
   }' |\
cat >! phenix_opts_fixangle.txt

egrep "^extern" refmac_opts_fixtorsion.txt |\
awk '$2=="torsion"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     c3=$23;r3=$25;a3=$27;f3=$29;\
     c4=$32;r4=$34;a4=$36;f4=$38;\
     v=$40;s=$42;\
   print "    #",c1,c2,c3,c4,a1,a2,a3,a4,r1,r2,r3,r4,f1,f2,f3,f4,v,s;\
   print "    dihedral {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      atom_selection_3 = \"name",a3,"and resseq",r3,"and chain",c3,"and altid",f3 "\"";\
   print "      atom_selection_4 = \"name",a4,"and resseq",r4,"and chain",c4,"and altid",f4 "\"";\
   print "      angle_ideal =",v;\
   print "      sigma =",s;\
   print "    }";\
   }' |\
cat >! phenix_opts_fixtorsion.txt

cat phenix_opts_*.txt |\
awk 'BEGIN{print "  geometry_restraints.edits {"} {print} END{print "  }"}' |\
awk 'BEGIN{print "refinement {"} {print} END{print "}"}' |\
cat >! phenix_opts_add.eff


skipfudge:









# keep record of worst of each type
echo -n "" >! ${t}_badatoms.txt
foreach type ( bond angle torsion plane chiral nonbond cbetadev )
  set thresh = 2.5
  awk -v thresh=2.5 '$2>thresh{$1="";print}' ${t}_worst${type}.txt |\
  awk -F "|" 'NF==2{print $2,$1+0}' |\
  awk -v type=$type '{for(i=1;i<NF;++i)print $i,$NF,type}' |\
  tee -a ${t}_badatoms.txt > /dev/null
end
sort -k2gr ${t}_badatoms.txt >! ${outprefix}_worstatoms.txt
awk '{print $1}' ${outprefix}_worstatoms.txt >! worstatoms_fliplist.txt

# one-liner of Mp stats
cat ${outprefix}_molprobity.log ${t}_worsts.txt |\
awk '/MolProbity score/{Mp=$NF}\
   /Clashscore   /{Cl=$NF}\
   $NF=="bond"{wbn=$1}\
   $NF=="angle"{wag=$1}\
   $NF=="torsion"{wdh=$1}\
   $NF=="plane"{wpl=$1}\
   $NF=="chiral"{wch=$1}\
   $NF=="nonbond"{nb_e=$1}\
   $NF=="omega"{om=sqrt(($(NF-1)-180)^2)}\
   $NF=="rama"{ra=$1}\
   $NF=="rota"{ro=$1}\
   $NF=="CBdev"{CB=$1}\
   $NF=="clash"{cle=$1}\
   END{print Mp+0,Cl+0,nb_e,wbn,wag,wdh,wch,wpl,om,ra,ro,CB;}' |\
cat >! ${t}.txt
set MPstats = `cat ${t}.txt`

set worstats = ( $worst_omega[1] $worst_rama[1] $worst_rota[1]  $worst_CB[1] )
set avgstats = ( $avg_nonbond[1] $avg_bond[1] $avg_angle[1] $avg_torsion[1] $avg_chiral[1] $avg_plane[1] $avg_omega[1] $avg_rama[1] $avg_rota[1] $avg_CB[1] )

echo "$MPstats     $avgstats   $outprefix" |\
tee ${outprefix}_stats.txt 

# Mp Cl nb_E wbond wang wdihe wchir wlpane     om  ra ro CB  nonbond bond angle tor chir plane omega rama rot CB pdb
#  1 2   3    4     5    6     7      8         9  10 11 12  13      14   15    16  17   18    19    20   21  22 23




# now try to incorporate R factor?
if(! -s  ${t}Rstats.txt) then
  echo "0 n/d n/d n/d n/d n/d n/d n/d n/d n/d n/d n/d 0 R" >! ${t}Rstats.txt
  # trial Rw Rf FOM LL LLf rmsBond Zbond rmsAngle zAngle rmsChiral function vdw   trial    "R"
  #   1   2  3  4   5  6   7       8     9        10     11         12       13   14       15
endif
if(! -s ${t}movement.txt) then
  echo "x moved x A x occ x B x x d" >! ${t}movement.txt
endif


awk '{print $0,"M"}' ${outprefix}_stats.txt >! ${t}Mstats.txt
# Mp Cl  nb_E wbond wang wdihe wchir wlpane     om  ra ro CB  nonbond bond angle tor chir plane omega rama rot CB pdb
#  1 2   3    4     5    6     7      8          9  10 11 12  13      14   15    16  17   18    19    20   21  22 23

echo $energy |\
cat - ${t}movement.txt ${t}Rstats.txt ${t}Mstats.txt |\
awk 'NR==1{energy=$1;next}\
    $NF=="d"{dxyz=$3;dB=$7;next}\
    $NF=="R"{Rw=$2;Rf=$3;bond=$7;ang=$9;vdw=$13;\
       rE=((Rf-2.0)/2.0)^2;\
       rstats=$2" "$3" "$7" "$9" "$11" "$13;next}\
    $NF=="M"{Mp=$1;Cl=$2; score=rE+energy;\
    $NF="";\
    if(rstats)print score,energy,rstats,$0,dxyz,dB,fliplist}' |\
awk 'NF>=31' >! ${outprefix}_score.txt
# score energy Rw Rf bond ang chir vdw   Mp Cl  nb_e wbond wang wdihe wchir wlpane om  ra ro CB  avnb  avbond avangle avtor avchir avplane avomega avrama avrot avCB  outprefix  dxyz dB  
  # 1     2      3  4  5    6   7    8     9  10  11   12    13   14    15    16     17  18 19 20  21    22     23      24    25     26      27      28     29    30    31         32   33

cat ${outprefix}_score.txt


exit:
if($?BAD) then
    echo "ERROR: $BAD"
    exit 9
endif

#cp ${t}.geo ${outprefix}.geo
#cp ${t}geo.txt ${outprefix}noH.geo


if( ! $debug && "$tempfile" != "" && "$tempfile" != "./") then
   rm -rf ${tmpdir}/
   rm -f ${t}*
endif

exit




