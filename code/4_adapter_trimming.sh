#!/bin/bash

# Illumina universal adapters revealed in fastqc

WD=$( pwd )
FQDIR=${WD}/data/fqfiles/
SAMPLELIST=( ${FQDIR}/*/ )
SCRIPTSDIR=${WD}/code/cutadapt/
OUTDIR=${WD}/output/cutadapt/

# make dirs
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR}

# loop
for f in "${SAMPLELIST[@]}" ; do

SAMPLE=$( basename "${f}" )

RUNSCRIPT=${SCRIPTSDIR}/qsub/cutadapt.${SAMPLE}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load python/2.7.13

cutadapt -m 20 --max-n 0.15 -a AGATCGGAAGAG -A AGATCGGAAGAG -o ${OUTDIR}/${SAMPLE}_1.trimmed.fq.gz -p ${OUTDIR}/${SAMPLE}_2.trimmed.fq.gz ${FQDIR}/${SAMPLE}/${SAMPLE}_*_1.fq.gz ${FQDIR}/${SAMPLE}/${SAMPLE}_*_2.fq.gz

EOF

qsub -l ncpus=4,mem=4G,walltime=3:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N trim_${SAMPLE} ${RUNSCRIPT}

done

