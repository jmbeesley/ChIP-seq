#!/bin/bash

# input sample
# scale samples

WD=$( pwd )
INPUTDIR=${WD}/output/filtering/filteredbams/
SCRIPTSDIR=${WD}/code/signal_tracks/
OUTDIR=${WD}/output/signal_tracks/
SAMPLES=( ${INPUTDIR}/*.nodup.sorted.bam )
CONTROL=$( grep control ${WD}/docs/sample_info.txt | cut -f1 )
THREADS=8   ## check

mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR} 

# loop
for f in "${SAMPLES[@]}" ; do

SAMPLE=$( basename "${f%.nodup.sorted.bam}" )
PREFIX=${SAMPLE}
FINALBAM=${INPUTDIR}/${PREFIX}.nodup.sorted.bam  ### final bam file
FINALBW=${OUTDIR}/${PREFIX}.coverage.bw 
CONTROLBAM=${INPUTDIR}/${CONTROL}.nodup.sorted.bam
LOG2RATIOBW=${OUTDIR}/${PREFIX}.l2rcontrol.bw 

RUNSCRIPT=${SCRIPTSDIR}/qsub/maketracks.${PREFIX}.sh

cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load python/2.7.13

# Create coverage track
bamCoverage -b ${FINALBAM} \
--normalizeUsing RPGC \
--effectiveGenomeSize 2864785220 \
--binSize 10 \
--extendReads 200 \
-o ${FINALBW}

# compare to input
bamCompare -b1 ${FINALBAM} \
-b2 ${CONTROLBAM} \
-p ${THREADS} \
--extendReads 200 \
-bs 10 \
-o ${LOG2RATIOBW}

EOF

qsub -l ncpus=10,mem=48G,walltime=12:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N maketracks_${PREFIX} ${RUNSCRIPT}

done


