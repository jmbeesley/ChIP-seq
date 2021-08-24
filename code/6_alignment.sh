
## run on hpcpbs01
## Histone mark ChIP-seq
## ENCODE ChIP-seq alignment parameters for bowtie2

DATADIR=output/cutadapt/
FQFILES=( ${DATADIR}/*_1.trimmed.fq.gz )
BT2IDX=/reference/genomes/GRCh37_ICGC_standard_v2/indexes/BOWTIE2_2.2.9/GRCh37_ICGC_standard_v2.fa
SCRIPTSDIR=code/alignment/
OUTDIR=output/alignment/bowtie2/
BAMSTATSDIR=output/alignment/bowtie2stats/
THREADS=8 ## check

# make dirs
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR} ${BAMSTATSDIR}

# loop
for f in "${FQFILES[@]}" ; do

SAMPLE=$( basename "${f%_1.trimmed.fq.gz}" )
PREFIX=${SAMPLE}
FASTQ1=${DATADIR}/${PREFIX}_1.trimmed.fq.gz
FASTQ2=${DATADIR}/${PREFIX}_2.trimmed.fq.gz
OUTPUTBAM=${OUTDIR}/${PREFIX}.bam
BAMINDXSTATS=${BAMSTATSDIR}/${PREFIX}.idxstats
BAMFLAGSTATS=${BAMSTATSDIR}/${PREFIX}.flagstat
LOG=${SCRIPTSDIR}/log/${PREFIX}.log
RUNSCRIPT=${SCRIPTSDIR}/qsub/bt2_alignment.${PREFIX}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load bowtie2/2.2.9
module load samtools/1.9

# align reads
bowtie2 --threads ${THREADS} -x ${BT2IDX} --maxins 2000 -1 ${FASTQ1} -2 ${FASTQ2} 2> ${LOG} |
        samtools view -Su /dev/stdin |
        samtools sort -o ${OUTPUTBAM} -

# index bam file
samtools index ${OUTPUTBAM}

# stats
samtools idxstats ${OUTPUTBAM} > ${BAMINDXSTATS}
samtools flagstat -@ ${THREADS} ${OUTPUTBAM} > ${BAMFLAGSTATS}

EOF

qsub -l ncpus=16,mem=24G,walltime=24:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N align_${SAMPLE} ${RUNSCRIPT}

done
