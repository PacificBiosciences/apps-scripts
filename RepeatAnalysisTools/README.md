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
[sample].[label].repeatCount_kde.png

![RepeatCount plot Example](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/bc1002.FMR1.repeatCount_kde.png)

[sample].[label].waterfall.png

![Waterfall plot Example](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/bc1002.FMR1.waterfall.png)

[sample].[label].repeatCounts.csv

    $ column -ts, bc1002.FMR1.repeatCounts.csv | head
    readName                           AGG  CGG
    m54006_180727_184845/10027338/ccs  1    200
    m54006_180727_184845/10486148/ccs  1    141
    m54006_180727_184845/10945324/ccs  1    245
    m54006_180727_184845/12452147/ccs  1    23
    m54006_180727_184845/12649105/ccs  1    22
    m54006_180727_184845/12649463/ccs  2    132
    m54006_180727_184845/13959293/ccs  1    213
    m54006_180727_184845/15335794/ccs  1    23
    m54006_180727_184845/15336396/ccs  1    22

[sample].[label].summary.csv

    $ column -ts, bc1002.FMR1.summary.csv
    totalReads     2240
    spanningReads  165
    oneSided       84
    poorAlignment  1991
    reference      rources/FMR1_L446_R503.fasta

## minimap2_e40.sh
Map repeat-extension data to targets, parameterized to reduce alignment cost for long expansions

### Dependencies
 - minimap2
 - samtools

### Usage
    minimap2_e40.sh REFFASTA QUERY OUTBAM
