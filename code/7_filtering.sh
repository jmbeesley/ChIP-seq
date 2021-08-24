#!/bin/bash

# Post-alignment Filtering of uninformative reads: PCR duplicates (picard), non-uniuqe alignments,

WD=$( pwd )
INPUTDIR=${WD}/output/alignment/bowtie2/
SCRIPTSDIR=${WD}/code/filtering/
OUTDIR=${WD}/output/filtering/filteredbams/
SAMPLES=( ${INPUTDIR}/*.bam )
BAMSTATSDIR=${WD}/output/filtering/filterstats/
THREADS=8   ## check

mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR} ${BAMSTATSDIR}

# loop
for f in "${SAMPLES[@]}" ; do

SAMPLE=$( basename "${f%.bam}" )
PREFIX=${SAMPLE}
RAWBAMFILE=${INPUTDIR}/${PREFIX}.bam
NSORTEDRAWBAMFILE=${OUTDIR}/${PREFIX}.nsorted.bam
TEMPBAM=${OUTDIR}/${PREFIX}.temp.bam
FIXMATEBAM=${OUTDIR}/${PREFIX}.fixmate.bam
DUPMARKBAM=${OUTDIR}/${PREFIX}.dupmark.bam
DUPREMBAM=${OUTDIR}/${PREFIX}.nodup.bam
FINALBAM=${OUTDIR}/${PREFIX}.nodup.sorted.bam  ### final bam file
FINALBW=${OUTDIR}/${PREFIX}.coverage.bw  ### coverage bigwig file
FINAL_NMSRT_BAM=${OUTDIR}/${PREFIX}.nodup.nsorted.bam
FINAL_BEDPE_FILE=${OUTDIR}/${PREFIX}.nsorted.bedpe.gz ### BEDPE file (with read pairs on each line) 
FINAL_TA_FILE=${OUTDIR}/${PREFIX}.SE.tagAlign.gz  ### tagAlign file (virtual single end)
BAMIDXSTATS=${BAMSTATSDIR}/${PREFIX}.dupmark.idxstats ### for chrM reads
DUPMETRICS=${BAMSTATSDIR}/${PREFIX}.dupmetrics.log
DUPMARKFLAGSTATS=${BAMSTATSDIR}/${PREFIX}.dupmark.flagstat
FINALBAMIDXSTATS=${BAMSTATSDIR}/${PREFIX}.nodup.sorted.idxstats
FINALBAMFLAGSTATS=${BAMSTATSDIR}/${PREFIX}.nodup.sorted.flagstat
RUNSCRIPT=${SCRIPTSDIR}/qsub/filter.${PREFIX}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load picard/2.19.0
module load samtools/1.9
module load bedtools/2.29.0
module load python/2.7.13

# mapping quality > 30, unmapped reads, mate unmapped, not primary alignment, reads failing platform and unplaced contigs. Only keep properly paired reads
# Will produce name sorted BAM
samtools view  -@ ${THREADS} -F 1804 -f 2 -q 30 -b ${RAWBAMFILE} |
samtools sort -n - -o ${NSORTEDRAWBAMFILE}

# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools fixmate -r ${NSORTEDRAWBAMFILE} ${TEMPBAM}
samtools view -F 1804 -f 2 -u ${TEMPBAM} | samtools sort - -o ${FIXMATEBAM}
rm -f ${TEMPBAM}

# Mark duplicates with Picard
java -XX:ParallelGCThreads=8 -jar /software/picard/picard-tools-2.19.0/picard.jar MarkDuplicates INPUT=${FIXMATEBAM} OUTPUT=${DUPMARKBAM} METRICS_FILE=${DUPMETRICS} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

# Index
samtools index ${DUPMARKBAM}

# # idxstats and flagstat for BAM with duplicates
# samtools idxstats ${DUPMARKBAM} > ${BAMIDXSTATS}
# samtools flagstat -@ ${THREADS} ${DUPMARKBAM} > ${DUPMARKFLAGSTATS}

# Filter duplicates. Index final position sorted BAM. Create final name sorted BAM
samtools view -@ ${THREADS} -F 1804 -f 2 -q 30 -b ${DUPMARKBAM} > ${DUPREMBAM}
samtools sort -@ ${THREADS} -o ${FINALBAM} ${DUPREMBAM}

# Index Final BAM file
samtools index ${FINALBAM}
samtools idxstats ${FINALBAM} > ${FINALBAMIDXSTATS}
samtools flagstat -@ ${THREADS} ${FINALBAM} > ${FINALBAMFLAGSTATS}

# Sort by read name
samtools sort -n ${FINALBAM} -o ${FINAL_NMSRT_BAM}

# Create BEDPE file
bedtools bamtobed -bedpe -mate1 -i ${FINAL_NMSRT_BAM} | gzip -nc > ${FINAL_BEDPE_FILE}

# Create TagAlign file
zcat ${FINAL_BEDPE_FILE} | 
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' | 
gzip -nc > ${FINAL_TA_FILE}

EOF

qsub -l ncpus=10,mem=80G,walltime=12:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N filter_${SAMPLE} ${RUNSCRIPT}

done


