#! /bin/bash

OUTDIR=$1
JOBNOS="${@:2}"

if [ -z $1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then 
 echo -e "\nUsage: $(basename $0) outputDir jobnumber1 [jobnumber2 jobnumber3 ...]\n"
 exit 0
fi

for job in ${JOBNOS}; do
 #Sample name from SM value in RG
 sample=$(samtools view -H $(pbservice get-job --show-files ${job} \
                             | awk '/mapped.bam / {print $NF}') \
         | grep '@RG' | tr '\t' '\n' | awk -F: '/SM/ {print $2}')

 echo "Copying data from job ${job}, sample ${sample}"

 #Copy both bam and index
 pbservice get-job --show-files ${job} | awk '/mapped.bam/ {print $NF}' | parallel cp {} ${OUTDIR}/${sample}.{/}
done
