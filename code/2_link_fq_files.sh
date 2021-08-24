#!/bin/bash

WD=$( pwd )
DATADIR=data/raw_data/
SAMPLE=( ${DATADIR}/*/raw_data/*/ )
FQFILES=( ${DATADIR}/*/raw_data/*/*.fq.gz  )
FQDIR=data/fqfiles/

# create dir and put links into sample dir
for i in ${SAMPLE[@]} ; do

    sample=$( basename ${i}  )
    mkdir -p ${FQDIR}/${sample} 
    ( 
      cd ${FQDIR}/${sample}
      ln -s ${WD}/${DATADIR}/*/raw_data/${sample}/*.fq.gz . 
    )

done