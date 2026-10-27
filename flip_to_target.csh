#! /bin/tcsh -f
#
#  flip x-ray equivalent atoms so they match a target template
#
#  2026-10-27 (claude_numbers, JH: "make flip-to-target fix chi5"): added to this script: ARG is never renamed to match the
#  target; every ARG conformer is named canonically from its own geometry (chi5 = CD-NE-CZ-NH1 beyond +/-90 deg ->
#  NH1<->NH2, HH11<->HH21, HH12<->HH22; positions unchanged; logged "chi5 fix: ..." on stderr).  phenix.refine with
#  flip_symmetric_amino_acids=False cannot repair a non-canonical guanidinium (chi5 restraint ideal 0, period 1).
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

# ARG: canonical NH1/NH2 naming from each conformer's own chi5 (CD-NE-CZ-NH1); blank-altloc CD/NE/CZ serve every conformer
awk 'function dih(a,b,c,d,  b0x,b0y,b0z,b1x,b1y,b1z,b2x,b2y,b2z,n,vx,vy,vz,wx,wy,wz,t,cx,cy,cz){\
       b0x=a[1]-b[1];b0y=a[2]-b[2];b0z=a[3]-b[3]; b1x=c[1]-b[1];b1y=c[2]-b[2];b1z=c[3]-b[3]; b2x=d[1]-c[1];b2y=d[2]-c[2];b2z=d[3]-c[3];\
       n=sqrt(b1x*b1x+b1y*b1y+b1z*b1z); b1x/=n;b1y/=n;b1z/=n;\
       t=b0x*b1x+b0y*b1y+b0z*b1z; vx=b0x-t*b1x;vy=b0y-t*b1y;vz=b0z-t*b1z;\
       t=b2x*b1x+b2y*b1y+b2z*b1z; wx=b2x-t*b1x;wy=b2y-t*b1y;wz=b2z-t*b1z;\
       cx=b1y*vz-b1z*vy;cy=b1z*vx-b1x*vz;cz=b1x*vy-b1y*vx;\
       return atan2(cx*wx+cy*wy+cz*wz, vx*wx+vy*wy+vz*wz)*57.29577951}\
  FNR==1{++pass}\
  /^ATOM|^HETAT/{typ=substr($0,18,3);atom=substr($0,13,4);gsub(" ","",atom);alt=substr($0,17,1);ch=substr($0,22,1);res=substr($0,23,4)+0}\
  pass==1{if((/^ATOM|^HETAT/) && typ=="ARG" && (atom=="CD"||atom=="NE"||atom=="CZ"||atom=="NH1")){\
            k=ch" "res" "alt" "atom; X[k,1]=substr($0,31,8)+0;X[k,2]=substr($0,39,8)+0;X[k,3]=substr($0,47,8)+0; have[k]=1;\
            if(atom=="NH1") nh1[ch" "res" "alt]=1}\
          next}\
  pass==2 && done==0{done=1;\
     for(key in nh1){split(key,q," ");ch0=q[1];res0=q[2];a0=substr(key,length(ch0" "res0" ")+1,1);ok=1;\
       for(i=1;i<=4;i++){nm=(i==1?"CD":(i==2?"NE":(i==3?"CZ":"NH1")));kk=ch0" "res0" "a0" "nm;if(have[kk]==0)kk=ch0" "res0"   "nm;\
         if(have[kk]==0){ok=0;break} for(j=1;j<=3;j++)P[i,j]=X[kk,j]}\
       if(ok==0)continue;\
       for(j=1;j<=3;j++){p1[j]=P[1,j];p2[j]=P[2,j];p3[j]=P[3,j];p4[j]=P[4,j]}\
       c5=dih(p1,p2,p3,p4);\
       if(c5>90||c5<-90){fix[key]=1;printf("chi5 fix: ARG %s %d alt \"%s\" chi5 %.1f -> NH1/NH2 exchanged\n",ch0,res0,a0,c5) > "/dev/stderr"}}}\
  pass==2{if((/^ATOM|^HETAT|^ANISOU/) && typ=="ARG" && fix[ch" "res" "alt]){nm=substr($0,13,4);\
            if(nm==" NH1")nm=" NH2";else if(nm==" NH2")nm=" NH1";\
            else if(nm=="HH11")nm="HH21";else if(nm=="HH21")nm="HH11";\
            else if(nm=="HH12")nm="HH22";else if(nm=="HH22")nm="HH12";\
            $0=substr($0,1,12) nm substr($0,17)}\
          print}' ${tempfile}out ${tempfile}out |\
cat >! ${tempfile}out2
mv ${tempfile}out2 $outfile

rmsd.awk $pdbfile $outfile

rm -f ${tempfile}*

