REF=$1
QUERY=$2
OUTBAM=$3
SAMPLE=$4

if [ -z $SAMPLE ]; then
 SAMPLE=""
fi

minimap2 -a \
         -x map-pb \
         --eqx \
         -L \
         -O 5,56 \
         -E 4,0 \
         --lj-min-ratio 0.2 \
         -B 5 \
         -A 2 \
         -r 2000 \
         -z 400,50 \
         -Y \
         --secondary=no \
         -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
        "${REF}" "${QUERY}" | \
samtools view -bS - | samtools sort - > ${OUTBAM}
samtools index ${OUTBAM}
