#!/bin/bash

peak1=${1?first peak}
peak2=${2?second peak}
outdir=${3?output dir}

mkdir -p ${outdir}

bedtools intersect -a $peak1 -b $peak2 > ${outdir}/$(basename ${peak1%%.*}).common.peak
common=$(wc -l ${outdir}/$(basename ${peak1%%.*}).common.peak)
echo "$common"

bedtools subtract -A -a $peak1 -b $peak2 > ${outdir}/$(basename ${peak1}).only.bed
left=$(wc -l ${outdir}/$(basename ${peak1}).only.bed)

echo "$left"

bedtools subtract -A -b $peak1 -a $peak2 > ${outdir}/$(basename ${peak2}).only.bed
right=$(wc -l ${outdir}/$(basename ${peak2}).only.bed)

echo "$right"
