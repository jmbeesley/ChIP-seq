#!/bin/bash

# call peaks using MACS2 for each sample

WD=$( pwd )
INPUTDIR=${WD}/output/filtering/filteredbams/
SCRIPTSDIR=${WD}/code/peakcalling/
OUTDIR=${WD}/output/peakcalling/
SAMPLES=( ${INPUTDIR}/*.nodup.sorted.bam )
BLACKLIST=${WD}/docs/hg19.chrom.sizes
CHROMSIZES=${WD}/docs/hg19.encode.blacklist.merged.bed

# create directories
mkdir -p ${SCRIPTSDIR}/{log,qsub} ${OUTDIR}

# loop through samples
for i in ${SAMPLES[@]} ; do

SAMPLE=$( basename ${i%.nodup.sorted.bam} )
BAMFILE=${i}
RUNSCRIPT=${SCRIPTSDIR}/qsub/macs2.${SAMPLE}.sh

# run script
cat > ${RUNSCRIPT} << EOF
#!/bin/bash

module load python/2.7.13
module load ucsctools/20160223
module load bedtools/2.29.0

# call peaks
macs2 callpeak -t ${BAMFILE} -f AUTO -n ${SAMPLE} --outdir ${OUTDIR} -g hs -p 1e-2 --nomodel --shift 0 --extsize 147 --keep-dup all -B --SPMR

# create bigwig signal track (-B in callpeak writes the bdg file)
sort -k1,1 -k2,2n ${OUTDIR}/${SAMPLE}_treat_pileup.bdg > ${OUTDIR}/${SAMPLE}_treat_pileup.sorted.bdg
bedGraphToBigWig ${OUTDIR}/${SAMPLE}_treat_pileup.sorted.bdg ${CHROMSIZES} ${OUTDIR}/${SAMPLE}.nonorm.bigwig

# filter out peaks in blacklisted regions
bedtools intersect -v -a ${OUTDIR}/${SAMPLE}_peaks.narrowPeak -b ${BLACKLIST} |
  awk 'BEGIN{OFS="\t"} {if (\$5>1000) \$5=1000; print \$0}' |
  grep -P 'chr[0-9XY]+(?!_)' |
  gzip -nc > ${OUTDIR}/${SAMPLE}.narrowPeak.blfiltered.gz

EOF

qsub -l ncpus=1,mem=22G,walltime=3:00:00 -e ${SCRIPTSDIR}/log -o ${SCRIPTSDIR}/log -N macs2_${SAMPLE} ${RUNSCRIPT}

done


