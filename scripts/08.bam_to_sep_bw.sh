#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?please provied sample file}
threads=${2-8}
bs=${3-50}

ALIGN_DIR="analysis/2-read-align"
BW_DIR="analysis/4-normliazed-bw"

mkdir -p $BW_DIR

exec 0< $samples
cd $ALIGN_DIR

while read sample
do
    bams=$(ls ${sample}*flt.bam | tr '\n' ' ')
    if [ ! -f ../$(basename ${BW_DIR})/${sample}.tmp.bam ]; then
    echo "samtools merge -f -@ ${threads} ../$(basename ${BW_DIR})/${sample}.tmp.bam $bams  &&
          samtools index ../$(basename ${BW_DIR})/${sample}.tmp.bam" | bash
    fi
done

cd ../../
exec 0< $samples

cd ${BW_DIR}
while read sample
do
    ls
    bamCoverage -b ${sample}.tmp.bam -o ${sample}_fwd.bw -of bigwig  \
      --filterRNAstrand forward --binSize ${bs} --normalizeUsing CPM --effectiveGenomeSize 2864785220 \
      --extendReads 150 -p ${threads} 2> ../log/${sample}_fwd.log
    bamCoverage -b ${sample}.tmp.bam -o ${sample}_rvs.bw -of bigwig  \
      --filterRNAstrand reverse --binSize ${bs} --normalizeUsing CPM --effectiveGenomeSize 2864785220 \
      --extendReads 150 -p ${threads} 2> ../log/${sample}_rvs.log
    rm -f ${sample}.tmp.bam ${sample}.tmp.bam.bai
done
