# Example Walkthrough for No-Amp Data Prep and Repeat Analysis

## Recalling Adapters

### Requirements
 * bam2bam software suite
    Temporarily unavailable.  Contact [PacBio Appslab](jharting@pacificbiosciences.com)
 * Adapter sequences fasta
 * Raw subreads/scraps bam files

Adapter recalling can be done with a single command.  This takes as input the subreads/scraps bams, as well as a filte of the adapter sequences.  Read data can be passes as BAMs (both subreads/scraps) or a subreadset.xml that includes scraps in the dataset.

    bam2bam -j 20 \
            -o recalled/m54006_181206_225403.tetraloop_recall \
            --adpqc \
            --adapters adapters.tetraloop.fasta \
            m54006_181206_225403.subreadset.xml

# Call CCS Reads
We recommend calling CCS consensus on the subreads for further analysis.

    ccs -j 20 \
        recalled/m54006_181206_225403.tetraloop_recall.subreadset.xml \
        ccs/m54006_181206_225403.tetraloop_recall.consensusreadset.xml

## Alignment to Reference
We provide two bash scripts to parameterize alignment of CCS reads to the reference.  Both scripts use minimap2 as the aligner and have parameters set such that reads with extended repeat motifs are correctly mapped to the reference.

    pbmm2_extension.sh human_hs37d5.fasta ccs/m54006_181206_225403.tetraloop_recall.consensusreadset.xml \
                       aligned/m54006_181206_225403.tetraloop_recall.consensusalignmentset.xml \
                       [samplename]

## Restriction Diagnostics

No-Amp prep can be diagnosed for adapter configuration and restriction enzyme cut sites.

    python NoAmpRestrictionDiagnostics.py \
           -o diagnostics/ \
           -s recalled/m54006_181206_225403.tetraloop_recall.subreads.bam \
           -a adapters.tetraloop.fasta \
           aligned/m54006_181206_225403.tetraloop_recall.bam \
           human_hs37d5.targets_repeatonly.bed \
           human_hs37d5.fasta \
           restriction_enzymes.tsv

Inputs include a BED file of target regions (contig,start,stop,name,motifs):
    
    $ head resources/human_hs37d5.targets_repeatonly.bed
    4   3076554     3076717     HTT     CAG,CAA,CCG,CCA,CGG
    9   27573435    27573596    C9orf72 GGGGCC
    X   146993569   146993629   FMR1    CGG,AGG
    22  46191235    46191304    ATXN10  ATTCT,ATTCC,ATTTCT,ATTCCT

And a tsv file containing restriction sites to search for.  Degenerate motifs can be included by repeating restriction enzymes with different exact motifs.

    $ head resources/restriction_enzymes.tsv
    EcoR1   GAATTC
    BamH1   GGATCC
    EcoR5   GATATC
    Spe1    ACTAGT
    Avr2    CCTAGG
    BssSa1  CACGAG
    Acc1    GTAGAC
    Acc1    GTATAC
    Acc1    GTCGAC
    Acc1    GTCTAC

## On Target Counts

A quick count of on target molecules can be done with the recalled alignments file and a target BED.

    python countOnTarget.py \
           -o diagnostics/ \
           aligned/m54006_181206_225403.tetraloop_recall.bam \
           human_hs37d5.targets_repeatonly.bed

## Repeat Analysis

Basic repeat analysis to characterize the results can be accomplished with the RepeatAnalysisReport.py script.  This will generate visual and tabular results, as well it can be used to clip out and export the target regions from reads.

    python RepeatAnalysisReport.py \
           -o diagnostics/ \
           -s example \
           aligned/m54006_181206_225403.tetraloop_recall.bam \
           human_hs37d5.targets_repeatonly.bed \
           human_hs37d5.fasta
