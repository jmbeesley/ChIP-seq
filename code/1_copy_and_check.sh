#!/bin/bash

# md5 check
DATADIR=data/
CHECK=data/raw_data/md5_check/
SAMPLE=( data/raw_data/*/raw_data/*/ )

# create destination dir
mkdir -p ${CHECK}

# copy md5 txt files
for d in ${SAMPLE[@]} ; do

    sample=$( basename ${d} )
    echo checking "${sample}"
    
    md5sum ${d}*.fq.gz | sed "s@${d}@@" > ${CHECK}/${sample}_MD5.tmp
    
    # check match
    awk 'NR==FNR { a[$2]=$1 ; next }{ print $2"\t"$1"\t"a[$2] }'  ${d}/MD5.txt ${CHECK}/${sample}_MD5.tmp |
    awk '{ a = ($2==$3?"TRUE":"FALSE") ; print $0"\t"a }' > ${CHECK}/${sample}_MD5check.txt
  
    rm -f ${CHECK}/${sample}_MD5.tmp
done

echo "Done." 

cat ${CHECK}/*_MD5check.txt
