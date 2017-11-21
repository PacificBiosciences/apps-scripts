source /mnt/software/Modules/current/init/bash
module load parallel
module load samtools

USAGE="Usage: `basename $0` [SMRTLinkJobDir1] [SMRTLinkJobDir2] [nproc]"

if ! [[ $# == 3 ]]; then
        echo $USAGE
        exit 0
fi

if [ "$1" == "-h" ]; then
        echo $USAGE
        exit 0
fi


#SMRT Link Job dirs
# e.g. /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/userdata/jobs_root/076/076210
DIR1=$1
DIR2=$2
# number of procs
# e.g. number of pbalign jobs
nproc=$3

# file of filenames for mapped.alignmentset.bam 
find $DIR1 -name "mapped.alignmentset.bam" > BAM1.fofn
find $DIR2 -name "mapped.alignmentset.bam" > BAM2.fofn

# make/empty temp dirs
if [ -d "tmp1" ]; then
	rm tmp1/*
else
	mkdir tmp1
fi

if [ -d "tmp2" ]; then
        rm tmp2/*
else
        mkdir tmp2
fi


# get read IDs for each pbalign job
parallel --verbose -j 24 --files --tmpdir ./tmp1/ 'samtools view {} | cut -f1'  :::: BAM1.fofn
parallel --verbose -j 24 --files --tmpdir ./tmp2/ 'samtools view {} | cut -f1'  :::: BAM2.fofn


# concatenate read IDs, sort, uniq
> readIDs1.txt
for f in $(find tmp1/ -name "par*.par");
do
	cat $f >> readIDs1.txt
done

sort readIDs1.txt | uniq > tmp
mv tmp readIDs1.txt


> readIDs2.txt
for f in $(find tmp2/ -name "par*.par");
do
        cat $f >> readIDs2.txt
done

sort readIDs2.txt | uniq > tmp
mv tmp readIDs2.txt

comm -12 readIDs1.txt readIDs2.txt > sharedReads.txt
