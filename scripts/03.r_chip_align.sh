#!/bin/bash

set -e
set -u
set -o pipefail

# configuration
threads=8
index="index/hg19"
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
while read id
do
    if [ ! -f ${FQ_DIR}/$id.bt2.done ]
    then
    echo "bowtie2 --very-sensitive-local --mm -p $threads -x $index -U ${FQ_DIR}/$id.fastq.gz 2> ${LOG_DIR}/$id.bt2.log | \
    samtools sort -@ 2 -m 1G -T ${TMP_DIR}/${id} -o ${ALIGN_DIR}/${id}.sort.bam \
    && touch ${FQ_DIR}/$id.bt2.done 
    " | bash
    fi
done
