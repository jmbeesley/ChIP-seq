#!/bin/bash

# md5 check
CHECK=data/raw_data/md5_check/
SAMPLE=( data/raw_data/X*/raw_data/S*/ )

# create destination dir
mkdir -p ${CHECK}

# copy md5 txt files
for d in ${SAMPLE[@]}; do

    echo ${d}
    sample=$( basename ${d} )
    echo ${sample}
    ( 
      cd ${d}
      md5sum *.fq.gz | sed "s@${d}\/@@" > ${CHECK}/${sample}_MD5.tmp
    
      # check match
      awk 'NR==FNR { a[$2]=$1 ; next }{ print $2"\t"$1"\t"a[$2] }'  ${d}/MD5.txt ${CHECK}/${sample}_MD5.tmp |
      awk '{ a = ($2==$3?"TRUE":"FALSE") ; print $0"\t"a }' > ${CHECK}/${sample}_MD5check.txt
    
    )
    
      rm -f ${CHECK}/${sample}_MD5.tmp

done

# concatenate
cat ${CHECK}/*_MD5check.txt

