#!/bin/bash

DATADIR=data/raw_data/
FQFILES=( ${DATADIR}/X*/raw_data/S*/*.fq.gz  )
FQDIR=${WKDIR}/data/fqfiles/

# create dir and put links into sample dir
for i in ${FQFILES[@]} ; do

    mkdir -p ${FQDIR}/${i} 
    cd !$
    ln -s ${DATADIR}/X*/raw_data/${i}/*.fq.gz .

done