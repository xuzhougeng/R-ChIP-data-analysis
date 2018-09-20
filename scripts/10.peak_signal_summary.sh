#!/bin/bash

bed=${1?bed file}
fwd=${2?forwad bigwig}
rvs=${3?reverse bigwig}

exec 0< $bed

while read region
do
    chr=$( echo $region | cut -d ' ' -f 1)
    sta=$( echo $region | cut -d ' ' -f 2)
    end=$( echo $region | cut -d ' ' -f 3)
    (bigWigSummary $fwd $chr $sta $end 1 ; bigWigSummary $rvs $chr $sta $end 1 )| awk -v out=0 '{out=out+$1} END{print out/2}'
done 
