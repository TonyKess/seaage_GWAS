#!/bin/bash
#submit salmon jobs by chromosome 
while read -r x && read -r y <&3;
  do sed -i 's/'"$x"'/'"$y"'/g' $1 ;
  sbatch $1 ; 
  done <ssalist1 3<ssalist2
