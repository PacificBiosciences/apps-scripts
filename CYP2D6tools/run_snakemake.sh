#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

umask 002

BATCH=$1
BARCODES=$2
BIOSAMPLES=$3
HIFI=$4

if [ -z $3 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
 echo -e "\nUsage: $(basename $0) <batch_name> <barcode_fasta> <biosample_csv> <hifi_reads.bam>\n"
 exit 0
fi

mkdir -p "batches/${BATCH}/"
LOCKFILE="batches/${BATCH}/process_batch.lock"

# add lockfile to directory to prevent multiple simultaneous jobs
lockfile -r 0 "${LOCKFILE}" || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --configfile workflow/CYP2D6tools/config.yaml \
    --config batch="${BATCH}" \
             barcodeFa="${BARCODES}" \
             biosamples="${BIOSAMPLES}" \
             hifireads="${HIFI}" \
    --nolock \
    --local-cores 4 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend mamba \
    --latency-wait 120 \
    --cluster-config workflow/CYP2D6tools/cluster.yaml \
    --cluster "sbatch --partition={cluster.partition} \
                      --cpus-per-task={cluster.cpus} \
                      --output={cluster.out} {cluster.extra} " \
    --snakefile workflow/CYP2D6tools/Snakefile
