#! /bin/tcsh -f
#
#  flip x-ray equivalent atoms so they match a target template
#
#
set pdbfile = "$1"
set target = "$2"

if(! -e "$pdbfile") then
  cat << EOF
usage: $0 wrongflips.pdb target.pdb
EOF
  exit 9
endif

set outfile = flipped.pdb
set tempfile = tempfile_$$_

egrep "HIS|AS|GL|TYR|PHE" $pdbfile |\
awk '! /^ATOM|^HETAT/{next}\
 {atom=substr($0,12,5);gsub(" ","",atom)}\
 atom~/^[CON][DE]/{print}' |\
tee ${tempfile}flippable.pdb |\
awk '{typ=substr($0,18,3);p=0}\
  typ=="HIS"{++p;gsub(" CD2"," XXX");gsub(" ND1"," CD2");gsub(" XXX"," ND1");}\
  typ=="HIS"{++p;gsub(" CE1"," XXX");gsub(" NE2"," CE1");gsub(" XXX"," NE2");}\
  typ=="ASP"{++p;gsub(" OD1"," XXX");gsub(" OD2"," OD1");gsub(" XXX"," OD2");}\
  typ=="ASN"{++p;gsub(" OD1"," XXX");gsub(" ND2"," OD1");gsub(" XXX"," ND2");}\
  typ=="GLU"{++p;gsub(" OE1"," XXX");gsub(" OE2"," OE1");gsub(" XXX"," OE2");}\
  typ=="GLN"{++p;gsub(" OE1"," XXX");gsub(" NE2"," OE1");gsub(" XXX"," NE2");}\
  typ=="TYR" || typ=="PHE"{++p;gsub(" CD1"," XXX");gsub(" CD2"," CD1");gsub(" XXX"," CD2");}\
  typ=="TYR" || typ=="PHE"{++p;gsub(" CE1"," XXX");gsub(" CE2"," CE1");gsub(" XXX"," CE2");}\
  p{print}' |\
cat >! ${tempfile}flipped.pdb

rmsd.awk -v debug=1 ${tempfile}flippable.pdb $target |\
  awk '/moved/{c=substr($0,11,1);f=substr($0,6,1);r=substr($0,12,5);d=substr($0,25,9);\
     ssd[c" "f" "r]+=d*d}\
  END{for(cr in ssd)print cr,ssd[cr],"noflip"}' >! ${tempfile}noflip.txt
rmsd.awk -v debug=1 ${tempfile}flipped.pdb $target |\
awk '/moved/{c=substr($0,11,1);f=substr($0,6,1);r=substr($0,12,5);d=substr($0,25,9);\
     ssd[c" "f" "r]+=d*d}\
  END{for(cr in ssd)print cr,ssd[cr],"flip"}' >! ${tempfile}flip.txt

sort -k4g ${tempfile}flip.txt ${tempfile}noflip.txt |\
awk '{id=substr($0,1,9)} ! seen[id]{print;++seen[id]}' |\
awk '$NF=="flip"{print}' |\
cat - $pdbfile |\
awk '$NF=="flip"{c=substr($0,1,1);f=substr($0,3,1);r=substr($0,5,4)+0;\
    id=c" "f" "r;++flipme[id];next}\
! /^ATOM|^HETAT/{print;next}\
 {c=substr($0,22,1);f=substr($0,17,1);r=substr($0,23,4)+0;typ=substr($0,18,3);\
  id=c" "f" "r}\
 flipme[id]{;\
  if(typ=="HIS"){gsub(" CD2"," XXX");gsub(" ND1"," CD2");gsub(" XXX"," ND1");}\
  if(typ=="HIS"){gsub(" CE1"," XXX");gsub(" NE2"," CE1");gsub(" XXX"," NE2");}\
  if(typ=="ASP"){gsub(" OD1"," XXX");gsub(" OD2"," OD1");gsub(" XXX"," OD2");}\
  if(typ=="ASN"){gsub(" OD1"," XXX");gsub(" ND2"," OD1");gsub(" XXX"," ND2");}\
  if(typ=="GLU"){gsub(" OE1"," XXX");gsub(" OE2"," OE1");gsub(" XXX"," OE2");}\
  if(typ=="GLN"){gsub(" OE1"," XXX");gsub(" NE2"," OE1");gsub(" XXX"," NE2");}\
  if(typ=="TYR" || typ=="PHE"){gsub(" CD1"," XXX");gsub(" CD2"," CD1");gsub(" XXX"," CD2");}\
  if(typ=="TYR" || typ=="PHE"){gsub(" CE1"," XXX");gsub(" CE2"," CE1");gsub(" XXX"," CE2");}\
  if(typ=="TYR" || typ=="PHE"){gsub(" HD1"," XXX");gsub(" HD2"," HD1");gsub(" XXX"," HD2");}\
  if(typ=="TYR" || typ=="PHE"){gsub(" HE1"," XXX");gsub(" HE2"," HE1");gsub(" XXX"," HE2");}\
  atom=substr($0,12,5);gsub(" ","",atom);\
  Ee=substr(atom,1,1);$0=substr($0,1,77) Ee;\
 } {print}' |\
cat >! ${tempfile}out
mv ${tempfile}out $outfile

rmsd.awk $pdbfile $outfile

rm -f ${tempfile}*

