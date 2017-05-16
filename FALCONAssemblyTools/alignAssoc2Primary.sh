# Sarah B. Kingan
# 9 May 2017
# 
# Pacific Biosciences
# Applications Lab

USAGE="Usage: `basename $0` [FALCON_Assembly.fasta] [Primary Contig ID]"
EXAMPLE="Example: alignHaplotigs2Primary.sh myFCNAsm.fasta 000123F"

# if no arguments or help statement, print usage
if ! [[ $# == 2 ]]; then
	echo $USAGE
	echo $EXAMPLE
	exit 0
fi

if [ "$1" == "-h" ]; then
	echo $USAGE
	echo $EXAMPLE
	exit 0
fi


# reference file
ASM=$1
# primary contig name
PRIM_ID=$2

# is contigID correctly formatted?
REGEX="([0-9]{6}F)"
if ! [[ $PRIM_ID =~ $REGEX ]]; then
	echo "ERROR: incorrectly formatted contig ID"
	exit 0
fi

# check if samtools index file exists, create if doesn't
if ! [[ -s "$ASM.fai" ]]; then
	echo "samtools faidx $ASM"
	samtools faidx $ASM
fi
FAI=$ASM.fai

# get primary contig ID as formated in assembly file
P=$(grep $PRIM_ID $FAI | cut -f1 | grep -v '-')
HAPLOTIGS=$(grep $PRIM_ID $FAI | cut -f1 | grep '-')

# check is primary contigs exists
if [ -z "${P}" ]; then
	echo "NO PRIMARY CONTIG $P IN ASSEMBLY FILE $ASM"
	exit 0
fi

# check that there are haplotigs, otherwise, exit
if [ -z "${HAPLOTIGS}" ]; then
	echo "NO HAPLOTIGS FOR PRIMARY CONTIG $P"
	exit 0
fi


# fetch and format primary sequence
samtools faidx $ASM $P > ref.fa

# fetch and format haplotigs
for H in $HAPLOTIGS
do
	samtools faidx $ASM $H >> qry.fa
done

# run nucmer
nucmer -maxmatch ref.fa qry.fa -prefix $PRIM_ID
delta-filter -q $PRIM_ID.delta > $PRIM_ID.q.delta

