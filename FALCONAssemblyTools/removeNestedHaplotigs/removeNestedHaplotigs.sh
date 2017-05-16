# Sarah B. Kingan
# 7 April 2017
# 
# Pacific Biosciences
# Applications Lab

# if no arguments or help statement, print usage
if ! [ $1 ]; then
  echo "Usage: `basename $0` [FALCON-Unzip_Assembly.fasta]"
  exit 0
fi

if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` [FALCON-Unzip_Assembly.fasta]"
  exit 0
fi

# user input assembly file
ASM=$1

# check if samtools index file exists, create if doesn't
if ! [[ -s "$ASM.fai" ]]; then
    printf "samtools faidx $ASM\n"
    samtools faidx $ASM
fi
FAI=$ASM.fai

# list of primary contigs
PRIMARIES=$(grep -v '_' $FAI | cut -f1)

# to parse basename of contig and suffix (polishing annotation)
REGEX="([0-9]{6}F)(.*)"

# clear contig lists
> nested_haplotigs.txt
> retained_contigs_haplotigs.txt
> nested_haplotigs.fa
> retained_contigs_haplotigs.fa

# loop through primaries
for P in $PRIMARIES
do

	printf "\n\nPROCESSING PRIMARY CONTIG $P\n" >&2

	if [[ $P =~ $REGEX ]]
	then
		# parse contig name and suffix
		P_BASENAME=${BASH_REMATCH[1]}
		SUFFIX=${BASH_REMATCH[2]}

		echo $P_BASENAME >> retained_contigs_haplotigs.txt

		# list of halotigs associated with primary
		HAPLOTIGS=$(grep $P_BASENAME $FAI | cut -f1 | grep '_')


		if [ -z "${HAPLOTIGS}" ]; then
			printf "NO HAPLOTIGS FOR PRIMARY CONTIG $P\n" >&2


		else

			# make fasta file of primary for nucmer reference
			printf "GETTING PRIMARY CONTIG AS REFERENCE: samtools faidx $ASM $P"  >&2
			samtools faidx $ASM $P | sed 's/'$SUFFIX'//' > ref.fa

			> infile_list.txt

			for H in $HAPLOTIGS
			do

				printf "\nPROCESSING HAPLOTIG $H\n" >&2 

				#  make fasta file of haplotig
				H_BASENAME=$(echo $H | sed 's/'$SUFFIX'//')
				printf "GETTING HAPLOTIG AS QUERY: samtools faidx $ASM $H\n" >&2
				samtools faidx $ASM $H | sed 's/'$SUFFIX'//' > qry.fa

				# run nucmer
				printf "RUNNING NUCMER ALIGNMENT: nucmer -maxmatch ref.fa qry.fa -prefix $H_BASENAME\n" >&2
				nucmer -maxmatch ref.fa qry.fa -prefix $H_BASENAME

				# filter delta file
				printf "FILTERING ALIGNMENT FILE: delta-filter -g $H_BASENAME.delta > $H_BASENAME.global.delta\n" >&2
				delta-filter -g $H_BASENAME.delta > $H_BASENAME.global.delta

				# make coordinates file
				printf "GETTING ALIGNMENT COORDINATES: show-coords -qTl $H_BASENAME.global.delta > $H_BASENAME.coords\n" >&2
				show-coords -rTl $H_BASENAME.global.delta > $H_BASENAME.coords

				echo $H_BASENAME.coords >> infile_list.txt

				# clean up files
				rm $H_BASENAME.delta
				rm $H_BASENAME.global.delta

			done

			printf "\nFINDING NESTED HAPLOTIGS FOR $P\n" >&2

			# python script here with list of coords files
			./nestedHaplotigs.py infile_list.txt

			# clean up files
			for F in `cat infile_list.txt`
			do
				rm $F
			done


		fi
	fi
done


# clean up directory
rm ref.fa
rm qry.fa


# generate fasta files and indices
for i in `cat nested_haplotigs.txt`
do
	samtools faidx $ASM $i$SUFFIX >> nested_haplotigs.fa
done

for j in `cat retained_contigs_haplotigs.txt`
do
	samtools faidx $ASM $j$SUFFIX >> retained_contigs_haplotigs.fa
done

samtools faidx nested_haplotigs.fa
samtools faidx retained_contigs_haplotigs.fa

# print final summary
nNested=$(wc -l nested_haplotigs.txt | cut -f1 -d" ")
nRetained=$(wc -l retained_contigs_haplotigs.txt | cut -f1 -d" ")
lNested=$(cut -f2 nested_haplotigs.fa.fai |  awk '{ sum += $1 } END { print sum }')
lRetained=$(cut -f2 retained_contigs_haplotigs.fa.fai |  awk '{ sum += $1 } END { print sum }')

printf "****************************************************************\n"
printf "                 	FINAL SUMMARY                  		\n"
printf "  number of nested haplotigs removed: 		$nNested	\n"
printf "  length of sequence removed: 			$lNested	\n"
printf "  number of remaining contigs and haplotigs: 	$nRetained	\n"
printf "  length of remaining sequence: 		$lRetained	\n"
printf "\n"
printf "  REMOVAL COMPLETE.\n"
