source /mnt/software/Modules/current/init/bash
source /pbi/dept/bifx/gconcepcion/falcon/fc_env_170911/bin/activate
module load samtools

USAGE="Usage: `basename $0` [primaryContigsFasta] [haplotigsFasta]"

if ! [[ $# == 2 ]]; then
        echo $USAGE
        exit 0
fi

if [ "$1" == "-h" ]; then
        echo $USAGE
        exit 0
fi


# primary and haplotig fasta files
# supplied as command line args
P=$1
H=$2

# cleaned file names
PC=$(echo $P | sed 's/\.fa*/_cleaned\.fa/')
HC=$(echo $H | sed 's/\.fa*/_cleaned\.fa/')

# clean fasta (remove zero size contigs)
clean_fasta.py $P
clean_fasta.py $H

# index cleaned fasta
samtools faidx $PC
samtools faidx $HC

## split contigs into two files
## odd and even contig IDs in separate files
## haplotigs go with associated primary contig
# save contig IDs to file
cut -f1 $PC.fai > tmp_contigList.txt
cut -f1 $HC.fai >> tmp_contigList.txt
# base contig IDs (no haplotig or polish suffix)
cat tmp_contigList.txt | sed 's/_[0-9]\{3\}//g' | sort | uniq > tmp_baseContigList.txt
cat tmp_baseContigList.txt | sed 's/|arrow//g' > tmp
mv tmp tmp_baseContigList.txt
# split contig IDs into two files of similar size and length distribution
awk 'NR % 2 == 1' tmp_baseContigList.txt > tmp_baseContigList_even.txt
awk 'NR % 2 == 0' tmp_baseContigList.txt > tmp_baseContigList_odd.txt

# generate fasta files: contigset 1 (even line numbers)
echo "making contigSet1.fa"
> contigSet1.fa
for c in `cat tmp_baseContigList_even.txt`
do
        ID=$(grep $c $PC.fai | cut -f1)
        for i in $ID
        do
                samtools faidx $PC $i >> contigSet1.fa
        done
        ID=$(grep $c $HC.fai | cut -f1)
        for j in $ID
        do
                samtools faidx $HC $j >> contigSet1.fa
        done
done


# generate fasta files: contigset 2 (odd line numbers)
echo "making contigSet2.fa"
> contigSet2.fa 
for c in `cat tmp_baseContigList_odd.txt`
do
        ID=$(grep $c $PC.fai | cut -f1)
        for i in $ID
        do
                samtools faidx $PC $i >> contigSet2.fa
        done
        ID=$(grep $c $HC.fai | cut -f1)
        for j in $ID
        do   
                samtools faidx $HC $j >> contigSet2.fa
        done
done


