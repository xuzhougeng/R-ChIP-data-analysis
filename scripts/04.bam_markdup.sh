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
# mark duplication
while read id;
do
    if [ ! -f ${ALIGN_DIR}/${id}.mkdup.done ]
    then
    echo "sambamba markdup -t $threads ${ALIGN_DIR}/${id}.sort.bam ${ALIGN_DIR}/${id}.mkdup.bam \
    && touch ${ALIGN_DIR}/${id}.mkdup.done" | bash
    fi
done
