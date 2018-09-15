#!/bin/bash

sra_files=$1

sra_ids=$(egrep -o 'SRR[0-9]{5,9}' $sra_files)

mkdir -p analysis/log

for id in $sra_ids
do
    echo "prefetch $id >> analysis/log/download.log"
done
