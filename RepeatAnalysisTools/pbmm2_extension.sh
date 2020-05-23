REF=$1
QUERY=$2
OUTFILE=$3
SAMPLE=$4

nPROC=4

if [ -z $SAMPLE ]; then
 SAMPLE='""'
fi

pbmm2 align --sort \
            -j ${nPROC} \
            -r 10000 \
            -A 2 \
            -B 5 \
            -z 400 \
            -Z 50 \
            -o 5 \
            -O 56 \
            -e 4 \
            -E 0 \
            -L 0.1 \
            -c 10 \
            --sample ${SAMPLE} \
            $REF \
            $QUERY \
            $OUTFILE
