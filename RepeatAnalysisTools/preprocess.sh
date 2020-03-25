SUBREADS=$1
BARCODES=$2
REFERENCE=$3
OUTDIR=$4
NCHUNKS=$5
PROCPERCHUNK=$6
RESOURCE=$7

#Help/check inputs
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] || [ -z $5 ] || [ -z $6 ] || [ -z $7 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then 
 echo -e "\nUsage: $(basename $0) subreads.bam|xml barcodes.fasta|xml reference.fasta outdir nchunks nprocsPerChunk {local|cluster}"
 exit 0
fi

export REFERENCE
export PROCPERCHUNK
export SUBREADS
export OUTDIR
#check for chunking program
CHUNKCCS=$(dirname $0)/chunkCCS.sh
if [ ! -e $CHUNKCCS ]; then
    echo "The program chunkCCS.sh was not found in the same path as $(basename $0)"
    exit 1
fi

#check for executables
for prog in ccs lima pbmm2 parallel; do
    if [ ! $(type -P $prog) ]; then
        echo "Missing program $prog in command path"
        exit 1
    fi
done

main() {
    parallel mkdir -p $OUTDIR/{} ::: ccs demux align
    #run ccs
    $CHUNKCCS $NCHUNKS $PROCPERCHUNK $RESOURCE $(ccs_cmd)
    #run demux
    run_demux 
    #run align
    case $RESOURCE in
        local)
            cmd=align_cmd
            ;;
        cluster)
            cmd=cluster_cmd 
            ;;
        *)
            echo "Unknown resource.  {local|cluster}"
            exit 1
    esac       

    parallel --xapply $(eval "$cmd {1} {2} $(outbase align)") ::: $(outbase demux)*bam ::: $(ls $(outbase demux)*bam | parallel -k get_barcode)
}

outbase () {
    local subdir=$1
    local name=$(basename $SUBREADS | sed -e 's/subread/ccs/' -e 's/.bam$//' -e 's/.xml$//')
    echo "${OUTDIR}/$subdir/$name"
}

get_barcode() {
    local bam=$1
    basename $bam | awk -F. '{print $(NF-1)}'
}
export -f get_barcode

ccs_cmd () {
    echo ccs $SUBREADS $(outbase ccs).bam --disable-heuristics --draft-mode full
}

run_demux() {
    lima --ccs --same -j ${PROCPERCHUNK} \
         --split-bam-named --peek-guess  \
         $(outbase ccs).bam $BARCODES $(outbase demux).bam
}

align_cmd() {
    local qry=$1
    local smpl=$2
    local pref=$3
    echo "pbmm2 align --sort \
                -j ${PROCPERCHUNK} \
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
                -c 0 \
                --sample ${smpl} \
                $REFERENCE $qry ${pref}.${smpl}.bam"
}
export -f align_cmd

cluster_cmd() {
    local qry=$1
    local smpl=$2
    local pref=$3
    echo "qsub -cwd -pe smp ${PROCPERCHUNK} -b y -j y -V -S /bin/bash $(align_cmd $qry $smpl $pref)"
}

main
