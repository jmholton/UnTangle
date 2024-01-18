#! /bin/tcsh -f
#
#
#
cat << EOF >! README.md
# UnTangle Challenge
EOF

cat ~/html/challenge/twoconf/test.html |\
 awk '/img src/{gsub("img src=","img src=images/")} {print}' |\
  grep -v "TITLE" >! README.md 
  
git add . ; git commit -m "html stuff" ; git push

