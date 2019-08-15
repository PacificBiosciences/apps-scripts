# Example Walkthrough for No-Amp Data Prep and Repeat Analysis (manual)

## Demultiplex Subreads
The program *lima* is provided for demultiplexing.

    lima --same \
         --min-score 26 \
         --split-bam-named \
         --peek-guess \
         -j 20 \
         m54006_190802_093121.subreadset.xml \
         pacbio.barcodeset.xml \
         demux/m54006_190802_093121.subreadset.xml

## Call CCS Reads
We recommend calling CCS consensus and filtering to >=QV20 for further analysis.   

    parallel -j 2 ccs -j 12 --minPredictedAccuracy 0.99 {} ccs/{/.}.ccs.bam ::: demux/*bam

## Alignment to Reference
We provide two bash scripts to parameterize alignment of CCS reads to the reference.  Both scripts use minimap2 as the aligner and have parameters set such that reads with extended repeat motifs are correctly mapped to the reference.
    
    parallel -j 6 pbmm2_extension.sh human_hs37d5.fasta {} align/{/.}.aligned.bam {/.} ::: ccs/*bam

## Extract Target Region
Target regions are identified and clipped out by mapping sequence *flanking* the target coordinates used to each read and excising the sequence between.

    parallel python extractRegion.py {} \
                                     human_hs37d5.fasta \
                                     'X:146993569-146993629' \
                                     -o fastq/{/.}.extracted_FMR1.fastq ::: align/*bam 

## Generate Waterfall Plots
A useful and clear visualization.

    parallel python waterfall.py -m CGG,AGG -i {} -o reports/{/.}.waterfall.png ::: fastq/*fastq

## Generate Count Histograms
Another useful visualization.

    parallel python plotCounts.py -m CGG -i {} -o reports/{/.} ::: fastq/*fastq

## Generate Motif Count Table
Counts of exact string matches, as well as total length of sequence.

    parallel python countMotifs.py -m CGG,AGG -i {} -o reports/{/.}.counts.csv ::: fastq/*fastq
