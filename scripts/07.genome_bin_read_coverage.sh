#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?missing sample file}
chromsize=${2:-index/hg19.chrom.sizes}
size=${3:-3000}
ALIGN_DIR="analysis/2-read-align"
COV_DIR="analysis/3-genome-coverage"

mkdir -p ${COV_DIR}

exec 0< $samples

while read sample
do
    bedtools makewindows -g $chromsize -w $size | \
      bedtools intersect -b ${ALIGN_DIR}/${sample}.flt.bam -a - -c -bed > ${COV_DIR}/${sample}.ReadsCoverage
done
