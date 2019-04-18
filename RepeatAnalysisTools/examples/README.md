# Example Walkthrough for No-Amp Data Prep and Repeat Analysis (manual)

## Recalling Adapters

### Requirements
 * recalladapters 
 * Adapter sequences fasta
 * Raw subreads/scraps bam files

Adapter recalling can be done with a single command.  This takes as input the subreads/scraps bams, as well as a file of the adapter sequences.  Read data can be passes as BAMs (both subreads/scraps) or a subreadset.xml that includes scraps in the dataset.

    recalladapters -j 20 \
                   -o refarm/m54006_181206_225403.tetraloop_recall \
                   --adpqc \
                   --adapters adapters.tetraloop.fasta \
                   m54006_181206_225403.subreadset.xml

## Demultiplex Subreads
The program *lima* is provided for demultiplexing.  Note that SMRTbell templates in the No-Amp protocol have barcodes on only one side.

    lima --single-side \
         --min-score 26 \
         --split-bam-named \
         refarm/m54006_181206_225403.tetraloop_recall.subreadset.xml \
         barcodes_10.fasta \
         barcoded/m54006_181206_225403.barcoded.subredset.xml

## Call CCS Reads
We recommend calling CCS consensus and filtering to >=QV20 for further analysis.   

    ccs -j 20 \
        --minPredictedAccuracy 0.99 \
        barcoded/m54006_181206_225403.barcoded.BC1007--BC1007.subredset.xml \
        ccs/m54006_181206_225403.BC1007--BC1007.consensusreadset.xml

## Alignment to Reference
We provide two bash scripts to parameterize alignment of CCS reads to the reference.  Both scripts use minimap2 as the aligner and have parameters set such that reads with extended repeat motifs are correctly mapped to the reference.

    pbmm2_extension.sh human_hs37d5.fasta \
                       ccs/m54006_181206_225403.BC1007--BC1007.consensusreadset.xml \
                       aligned/m54006_181206_225403.BC1007--BC1007.consensusalignmentset.xml \
                       BC1007--BC1007 [optional]

## On Target Counts

A quick count of on target molecules can be done with the recalled alignments file and a target BED.

    python countOnTarget.py \
           -o diagnostics/ \
           aligned/m54006_181206_225403.tetraloop_recall.bam \
           human_hs37d5.targets_repeatonly.bed

## Repeat Analysis

Basic repeat analysis to characterize the results can be accomplished with the RepeatAnalysisReport.py script.  This will generate visual and tabular results, as well it can be used to clip out and export the target regions from reads.

    python RepeatAnalysisReport.py \
           -o report/ \
           -s BC1007--BC1007 \
           aligned/aligned/m54006_181206_225403.BC1007--BC1007.bam \
           human_hs37d5.targets_repeatonly.bed \
           human_hs37d5.fasta

Inputs include a BED file of target regions. *Columns* contig, start, stop, name, motifs. Motifs should be ordered such that the primary expected motif is first, and any potential interruption repeats follow.  Input BED files can contain a single target region.
    
    $ head resources/human_hs37d5.targets_repeatonly.bed
    4   3076554     3076717     HTT     CAG,CAA,CCG,CCA,CGG
    9   27573435    27573596    C9orf72 GGGGCC
    X   146993569   146993629   FMR1    CGG,AGG
    22  46191235    46191304    ATXN10  ATTCT,ATTCC,ATTTCT,ATTCCT
