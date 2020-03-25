#! /usr/bin/bash

DATADIR=$1
BED=$2
REFERENCE=$3

REPORTSCRIPT="$(dirname $(readlink -f "$0"))/RepeatAnalysisReport.py"

FLANKSIZE=100

for BAM in ${DATADIR}/align/*bam; do
 barcode="$(echo $BAM | awk -F. '{print $(NF-1)}')"
 outDir="$(readlink -f ${DATADIR})/reports/${barcode}"
 mkdir -p "${outDir}"
 echo ${barcode}
 python "${REPORTSCRIPT}" \
        -o "${outDir}" \
        -s "${barcode}" \
        -f ${FLANKSIZE} \
        "$(readlink -f ${BAM})" \
        "$(readlink -f ${BED})" \
        "$(readlink -f ${REFERENCE})"
done
