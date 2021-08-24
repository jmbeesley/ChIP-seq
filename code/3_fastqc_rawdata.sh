## run on job server

DATADIR=data/raw_data/
FQFILES=( ${DATADIR}/X*/raw_data/S*/*.fq.gz  )
SCRIPTSDIR=code/fastqc_raw/
OUTDIR=output/fastqc_raw/

# make dirs
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR}

# loop
for f in "${FQFILES[@]}" ; do

FILE=$( basename "${f%.fq.gz}" )
RUNSCRIPT=${SCRIPTSDIR}/qsub/fastqc_rawdata.${FILE}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load fastqc/0.11.5

# fastqc
fastqc -o ${OUTDIR} ${f}

EOF

qsub -l ncpus=1,mem=6G,walltime=3:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N fastqc1_${FILE} ${RUNSCRIPT}

done
