# Description
This repo contains miscellaneous stand-alone scripts for diagnosing and plotting NoAmp/Repeat Analysis data in CCS format.

UNDER CONSTRUCTION

## ccsRepeatAnalysisReport.py
Generate "waterfall" and repeat count kde plots for CCS reads.  
### Dependencies
 - matplotlib
 - numpy
 - pandas
 - seaborn
 - pbcore
 - mappy
### Usage
    $ python ccsRepeatAnalysisReport.py -h
    usage: ccsRepeatAnalysisReport.py [-h] [-t,--target TARGET]
                                      [-p,--preset PRESET] [-m,--motifs MOTIFS]
                                      [-o,--outDir OUTDIR] [-s,--sample SAMPLE]
                                      [-l,--label LABEL]
                                      ccsFastx
    
    generate repeat-kde plot,waterfall plot, and summary for repeat expansion ccs
    reads
    
    positional arguments:
      ccsFastx            input ccs reads, can be fasta or fastq
    
    optional arguments:
      -h, --help          show this help message and exit
      -t,--target TARGET  fasta reference of sequence flanking repeat region (2
                          sequences). Default resources/FMR1_L446_R503.fasta
      -p,--preset PRESET  preset motifs to search for. default 'FMR1'. Available:
                          FMR1,FUCHS,HTT,ALS
      -m,--motifs MOTIFS  Comma separated motifs to search for. Will over-ride
                          presets. Default none.
      -o,--outDir OUTDIR  Output directory. Default cwd
      -s,--sample SAMPLE  Sample name, prepended to output files. Default None
      -l,--label LABEL    Label, prepended to output files after sample name, e.g.
                          'FMR1'. Defaults to preset name, except if -m set

### Output
 - [sample].[label].repeatCount_kde.png
 - [sample].[label].waterfall.png
 - [sample].[label].repeatCounts.csv
 - [sample].[label].summary.csv

## minimap2_e40.sh
Map repeat-extension data to targets, parameterized to reduce alignment cost for long expansions

### Dependencies
 - minimap2
 - samtools

### Usage
    minimap2_e40.sh REFFASTA QUERY OUTBAM
