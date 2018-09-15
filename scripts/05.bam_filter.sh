#!/bin/bash

set -e
set -u
set -o pipefail

# configuration
threads=8
index=index/hg19
FQ_DIR="analysis/0-raw-data"
ALIGN_DIR="analysis/2-read-align"
LOG_DIR="analysis/log"
TMP_DIR="analysis/tmp"

mkdir -p ${ALIGN_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${TMP_DIR}

samples=$1

exec 0< $samples
# alignment
while read id;
do
    echo "samtools view -@ threads -bF 1804 -q 30 ${ALIGN_DIR}/${id}.mkdup.bam -o ${ALIGN_DIR}/${id}.flt.bam" | \
    qsub -V -cwd -pe openmpi $threads -N ${id}_mkdup -q wangjw -S /bin/bash -o /tmp -e /tmp
done
