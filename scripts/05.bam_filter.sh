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

samples=${1?missing sample files}

exec 0< $samples

# filteration
while read id;
do
    if [ ! -f ${ALIGN_DIR}/${id}.flt.done ]
    then
    echo "
    samtools view -@ threads -bF 1804 -q 30 ${ALIGN_DIR}/${id}.mkdup.bam -o ${ALIGN_DIR}/${id}.flt.bam\
    && samtools index ${ALIGN_DIR}/${id}.flt.bam \
    && touch ${ALIGN_DIR}/${id}.flt.done "| bash
    fi
done
