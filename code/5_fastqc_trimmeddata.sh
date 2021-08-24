#!/bin/bash


WD=$( pwd )
DATADIR=${WD}/output/cutadapt/
FQFILES=( ${DATADIR}/*.fq.gz  )
SCRIPTSDIR=${WD}/code/fastqc_trimmed/
OUTDIR=${WD}/output/fastqc_trimmed/

# make dirs
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR}

# loop
for f in "${FQFILES[@]}" ; do

FILE=$( basename "${f%.fq.gz}" )
RUNSCRIPT=${SCRIPTSDIR}/qsub/fastqc_trimmeddata.${FILE}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load fastqc/0.11.5

# fastqc
fastqc -o ${OUTDIR} ${f}

EOF

qsub -l ncpus=1,mem=6G,walltime=3:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N fastqc2_${FILE} ${RUNSCRIPT}

done
