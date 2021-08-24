#!/bin/bash

# fingerprint as a measure of ChIP enrichment

WD=$( pwd )
INPUTDIR=${WD}/output/filtering/filteredbams/
SCRIPTSDIR=${WD}/code/initial_qc/
OUTDIR=${WD}/output/initial_qc/
THREADS=1

mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR} 

# loop
for f in 1 ; do

BAMS=( ${INPUTDIR}/*.nodup.sorted.bam ) 
RUNSCRIPT=${SCRIPTSDIR}/qsub/initialqc.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load python/2.7.13

# plot fingerprints
plotFingerprint \
 -b ${BAMS[@]} \
--minMappingQuality 30 --skipZeros \
--numberOfSamples 50000 \
-T "Fingerprints of different samples"  \
--plotFile ${OUTDIR}/fingerprints.png \
--outRawCounts ${OUTDIR}/fingerprints.tab

EOF

qsub -l ncpus=1,mem=6G,walltime=2:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N initqc ${RUNSCRIPT}

done


