#! /bin/bash

TARGETBED=$1
REFERENCE=$2
OUTDIR=$3
BAMS="${@:4}"

if [ -z $1 ] || [ "$1" == "-h" ] || [ "$1" == "--help" ]; then 
 echo -e "\nUsage: $(basename $0) targets.BED reference.fasta outputDir aligned1.bam [aligned2.bam aligned3.bam ...]\n"
 exit 0
fi

#python instance to use
#change as needed
PYTHON=python

#scripts directory
GDIR=$(dirname $0)

#commands file
CMD="${OUTDIR}/runReportCmds.sh"

#fastq dir
FQDIR="${OUTDIR}/fastq"
#reports dir
RPDIR="${OUTDIR}/reports"
#cluster output dir
CLDIR="${OUTDIR}/cluster"

#make output dirs
parallel mkdir -p {} ::: "$FQDIR" "$RPDIR" "$CLDIR"

#make read counts
echo 'echo "Counting Reads"' > $CMD
echo "parallel ${PYTHON} ${GDIR}/countOnTarget.py {} ${TARGETBED} -s {/.} -o ${OUTDIR}/ ::: ${BAMS}" >> $CMD
echo -e "\n" >> $CMD

#make coverage plots
echo 'echo "Coverage Plots"' >> $CMD
echo "parallel ${PYTHON} ${GDIR}/coveragePlot.py {} -o ${OUTDIR}/{/.} -t ${TARGETBED} 1>/dev/null ::: ${BAMS}" >> $CMD
echo -e "\n" >> $CMD

while read -r chr start stop name motifs revcomp; do
 #use first motif for histograms
 pMotif=$(echo ${motifs} | cut -d, -f1)
 #revcomp flag
 if [ "$revcomp" == 1 ]; then
  rc="-r"
 else
  rc=""
 fi

 echo 'echo "'"${name} Extract Region\"" >> $CMD
 echo "parallel ${PYTHON} ${GDIR}/extractRegion.py ${rc} {} ${REFERENCE} \'${chr}:${start}-${stop}\' -o ${FQDIR}/{/.}.extracted_${name}.fastq ::: ${BAMS}" >> $CMD
 echo -e "\n" >> $CMD

 echo 'echo "'"${name} Waterfall Plots\"" >> $CMD
 echo "parallel ${PYTHON} ${GDIR}/waterfall.py -m ${motifs} -i {} -o ${RPDIR}/{/.}.waterfall.png 1>/dev/null ::: ${FQDIR}/*${name}*fastq" >> $CMD
 echo -e "\n" >> $CMD

 echo 'echo "'"${name} Histograms\"" >> $CMD
 echo "parallel ${PYTHON} ${GDIR}/plotCounts.py -m ${pMotif} -i {} -o ${RPDIR}/{/.} 1>/dev/null ::: ${FQDIR}/*${name}*fastq" >> $CMD
 echo -e "\n" >> $CMD

 echo 'echo "'"${name} Count Tables\"" >> $CMD
 echo "parallel ${PYTHON} ${GDIR}/countMotifs.py -m ${motifs} -i {} -o ${RPDIR}/{/.}.counts.csv ::: ${FQDIR}/*${name}*fastq" >> $CMD 
 echo -e "\n" >> $CMD

 echo 'echo "'"${name} Cluster Reads\"" >> $CMD
 echo "parallel ${PYTHON} ${GDIR}/clusterByRegion.py ${rc} -m ${motifs} -d -p ${CLDIR}/{/.}.${name} {} ${REFERENCE} \'${chr}:${start}-${stop}\' 1>/dev/null ::: ${BAMS}" >> $CMD
 echo -e "\n" >> $CMD
 
done < ${TARGETBED}

#Execute
/bin/bash $CMD
