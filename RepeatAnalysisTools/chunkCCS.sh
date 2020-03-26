   NCHUNKS=$1
PROCPERJOB=$2
  RESOURCE=$3       #'local' for local machine, 'cluster' to submit to cluster
    CCSCMD="${@:4}" #FORMAT: "ccs $inbam $outbam [options]"

MERGETHREADS=8

read ccs ibam obam opts <<< $CCSCMD
CCSprog=$(which $ccs)
BNAME=${obam%.*}

ccs_cmd() {
    local chunk=$1
    local total=$2
    local out=${BNAME}.chunk${chunk}.bam
    echo "${CCSprog} ${ibam} ${out} ${opts} --chunk ${chunk}/${total} -j ${PROCPERJOB}"
}

cluster_cmd() {
    local chunk=$1
    local total=$2
    echo "qsub -N chunkedCCS.chunk${chunk} -cwd -pe smp ${PROCPERJOB} -b y -j y -V -S /bin/bash $(ccs_cmd $chunk $total)"
}

get_count() {
    cnt=$(ls ${BNAME}.chunk*bam 2>>/dev/null | wc -l)
    if [ -z $cnt ]; then
        echo 0
    else
        echo $cnt
    fi
}

case $RESOURCE in
    local)
        cmd=ccs_cmd
        ;;
    cluster)
        cmd=cluster_cmd        
        ;;
    *)
        echo -n "Unknown resource; must be local or cluster"
        exit 0
esac

main() {
    echo "Starting jobs..."
    seq -w 1 ${NCHUNKS} | \
        parallel $(eval "$cmd {} $NCHUNKS")

    ndone=0
    while [ "$ndone" -lt "$NCHUNKS" ]; do
        newCount=$(get_count)
        if [ "$ndone" -lt "$newCount" ]; then
            ndone=$newCount
            echo "Finished chunks: $ndone"
        else
            sleep 5
        fi
    done
    
    samtools merge -@${MERGETHREADS} ${obam} ${BNAME}.chunk*bam
    if [ -e "$obam" ]; then
        rm ${BNAME}.chunk*bam*
        samtools index ${obam}
    fi

    echo "Done"
}
main 
