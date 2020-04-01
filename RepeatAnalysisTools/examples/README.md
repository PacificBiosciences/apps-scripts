# Example Walkthrough for No-Amp Data Prep and Repeat Analysis (manual)

## Call CCS Reads
We recommend calling CCS consensus with heuristics disabled, full draft mode, and filtering to >=QV20 (default) for further analysis.  Please note that disabling the heuristics produces more very long expanded reads but is considerably slower.  Chunking the raw subread data and processing on a cluster is recommended.

    $ ccs inSubreads.bam ccs/outCCS.bam -j 32 --disable-heuristics --draft-mode full

## Demultiplex Subreads
The program *lima* is provided for demultiplexing.

    lima --same                \
         --ccs                 \
         --split-bam-named     \
         --peek-guess          \
         -j 8                  \
         ccs/outCCS.bam        \
         pacbio.barcodeset.xml \
         demux/outDemuxed.bam

Note that bams will be split by barcode and barcode names will be infixed before the extension `.bam`.

## Alignment to Reference
We provide two bash scripts to parameterize alignment of CCS reads to the reference.  Both scripts use minimap2 as the aligner and have parameters set such that reads with extended repeat motifs are correctly mapped to the reference.
    
    parallel -j 6 pbmm2_extension.sh human_hs37d5.fasta {} align/{/.}.aligned.bam ::: demux/*bam

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

## Work directory
Primary contents of working directory (not including indice etc)

    tree -P '*bam|*fastq|*csv|*png'
    .
    ├── align
    │   ├── outDemuxed.bc1015--bc1015.aligned.bam
    │   ├── outDemuxed.bc1016--bc1016.aligned.bam
    │   ├── outDemuxed.bc1017--bc1017.aligned.bam
    │   ├── outDemuxed.bc1018--bc1018.aligned.bam
    │   └── outDemuxed.bc1019--bc1019.aligned.bam
    ├── ccs
    │   └── outCCS.bam
    ├── demux
    │   ├── outDemuxed.bc1015--bc1015.bam
    │   ├── outDemuxed.bc1016--bc1016.bam
    │   ├── outDemuxed.bc1017--bc1017.bam
    │   ├── outDemuxed.bc1018--bc1018.bam
    │   └── outDemuxed.bc1019--bc1019.bam
    ├── fastq
    │   ├── outDemuxed.bc1015--bc1015.aligned.extracted_FMR1.fastq
    │   ├── outDemuxed.bc1016--bc1016.aligned.extracted_FMR1.fastq
    │   ├── outDemuxed.bc1017--bc1017.aligned.extracted_FMR1.fastq
    │   ├── outDemuxed.bc1018--bc1018.aligned.extracted_FMR1.fastq
    │   └── outDemuxed.bc1019--bc1019.aligned.extracted_FMR1.fastq
    └── reports
        ├── outDemuxed.bc1015--bc1015.aligned.extracted_FMR1.counts.csv
        ├── outDemuxed.bc1015--bc1015.aligned.extracted_FMR1.insertSize.png
        ├── outDemuxed.bc1015--bc1015.aligned.extracted_FMR1.motifCount.png
        ├── outDemuxed.bc1015--bc1015.aligned.extracted_FMR1.waterfall.png
        ├── outDemuxed.bc1016--bc1016.aligned.extracted_FMR1.counts.csv
        ├── outDemuxed.bc1016--bc1016.aligned.extracted_FMR1.insertSize.png
        ├── outDemuxed.bc1016--bc1016.aligned.extracted_FMR1.motifCount.png
        ├── outDemuxed.bc1016--bc1016.aligned.extracted_FMR1.waterfall.png
        ├── outDemuxed.bc1017--bc1017.aligned.extracted_FMR1.counts.csv
        ├── outDemuxed.bc1017--bc1017.aligned.extracted_FMR1.insertSize.png
        ├── outDemuxed.bc1017--bc1017.aligned.extracted_FMR1.motifCount.png
        ├── outDemuxed.bc1017--bc1017.aligned.extracted_FMR1.waterfall.png
        ├── outDemuxed.bc1018--bc1018.aligned.extracted_FMR1.counts.csv
        ├── outDemuxed.bc1018--bc1018.aligned.extracted_FMR1.insertSize.png
        ├── outDemuxed.bc1018--bc1018.aligned.extracted_FMR1.motifCount.png
        ├── outDemuxed.bc1018--bc1018.aligned.extracted_FMR1.waterfall.png
        ├── outDemuxed.bc1019--bc1019.aligned.extracted_FMR1.counts.csv
        ├── outDemuxed.bc1019--bc1019.aligned.extracted_FMR1.insertSize.png
        ├── outDemuxed.bc1019--bc1019.aligned.extracted_FMR1.motifCount.png
        └── outDemuxed.bc1019--bc1019.aligned.extracted_FMR1.waterfall.png
