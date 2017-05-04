# creates two fasta files one of primary contigs, the other of haplotigs from FALCON-Unzip Assembly
# Written by Sarah Kingan

# script requires samtools
source /mnt/software/Modules/current/init/bash
module load samtools

# user input assembly file, must have ".fasta" file extension
ASM=$1

# check if samtools index file exists, create if doesn't
if ! [[ -s "$ASM.fai" ]]; then
    printf "samtools faidx $ASM\n"
    samtools faidx $ASM
fi
FAI=$ASM.fai

# get lists of primary and haplotigs
grep '_' $FAI | cut -f1 > h.txt
grep -v '_' $FAI | cut -f1 > p.txt

# assembly name without fasta ext
BASE=$(echo $ASM | sed 's/.fasta//')

# initialize files
> ${BASE}_primaryContigs.fasta
> ${BASE}_haplotigs.fasta

for p in `cat p.txt`
do 
        samtools faidx $ASM $p >> ${BASE}_primaryContigs.fasta
done

for h in `cat h.txt`
do 
        samtools faidx $ASM $h >> ${BASE}_haplotigs.fasta
done

# clean up
rm h.txt
rm p.txt
