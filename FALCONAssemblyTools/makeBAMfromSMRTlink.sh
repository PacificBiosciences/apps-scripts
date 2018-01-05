USAGE="Usage: `basename $0` [SMRTlinkJobTasksDir] [Reference] [Coords]"

# if no arguments or help statement, print usage
if ! [[ $# == 3 ]]; then
        echo $USAGE
        exit 0
fi

if [ "$1" == "-h" ]; then
        echo $USAGE
        exit 0
fi

# load samtools
source /mnt/software/Modules/current/init/bash
module load samtools
DIR=$1
REF=$2
COORDS=$3

# 24 chunks, alter as needed
for i in `seq 1 24`;
do
	samtools view -b -o tmp.$i.bam $DIR/pbalign.tasks.pbalign-$i/mapped.alignmentset.bam $COORDS
	samtools index tmp.$i.bam
done

files=$(find  . -name "*tmp*.bam" | sed 's/\.\///g' > bam.fofn)
samtools merge -b bam.fofn --reference $REF merged.bam
samtools index merged.bam 
samtools sort merged.bam > merged_sorted.bam
NAME=$(echo $COORDS | sed 's/:/_/')
samtools faidx $REF $COORDS > $NAME.fa
mv merged_sorted.bam $NAME.bam
samtools index $NAME.bam 
samtools faidx $NAME.fa

rm tmp*
rm merged*
rm bam.fofn
