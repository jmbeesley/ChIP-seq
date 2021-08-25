#!/bin/bash

# call peaks using MACS2 for each sample

WD=$( pwd )
INPUTDIR=${WD}/output/filtering/filteredbams/
SCRIPTSDIR=${WD}/code/peakcalling/
OUTDIR=${WD}/output/peakcalling/
SAMPLES=( ${INPUTDIR}/*.nodup.sorted.bam )

# new dirs
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR}

# loop thru samples
for i in ${SAMPLES[@]} ; do

SAMPLE=$( basename ${i%.nodup.sorted.bam} )
BAMFILE=${i}
RUNSCRIPT=${SCRIPTSDIR}/qsub/macs2.${SAMPLE}.sh

##
cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load python/2.7.13

#
macs2 callpeak -t ${BAMFILE} -f AUTO -n ${SAMPLE} --outdir ${OUTDIR} -g hs -p 1e-2 --nomodel --shift 0 --extsize 147 --keep-dup all -B --SPMR

EOF

qsub -l ncpus=1,mem=22G,walltime=3:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N macs2_${SAMPLE} ${RUNSCRIPT}

done