#!/bin/bash

file_in=$1

SRR_IDS=$( egrep -o 'SRR[0-9]{5,9}' ${file_in})
OUT_DIR="analysis/0-raw-data"
mkdir -p ${OUT_DIR}

for id in $SRR_IDS
do
   echo "fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@\$ac-\$si/\$ri' $id -O ${OUT_DIR} &"
done
