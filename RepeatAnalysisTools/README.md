# Description
This repo contains miscellaneous stand-alone scripts for diagnosing and plotting NoAmp/Repeat Analysis data in CCS format.

UNDER CONSTRUCTION

## countOnTarget.py
Generate table of ZMW counts per target.
### Dependencies
 - [pandas](https://pandas.pydata.org/)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
### Usage
    $ python countOnTarget.py -h
    usage: countOnTarget.py [-h] [-o,--outdir OUTDIR] inBAM inBED
    
    Generate table of ZMW counts per target
    
    positional arguments:
      inBAM               BAM file of aligned reads. Must have .bai index
      inBED               BED file with targets
    
    optional arguments:
      -h, --help          show this help message and exit
      -o,--outdir OUTDIR  directory to save output file. default cwd.

### Example
    $ python countOnTarget.py combined.consensusalignmentset.bam resources/human_hs37d5.targets.bed
    $ column -ts, onTargetCounts.tsv
    name     ctg  start      end        length  expected  onTargetZMWs  enrichment
    HTT      4    3076554    3076717    163     0.0172    1164          67683.5412
    C9orf72  9    27573435   27573596   161     0.0172    0             0.0000
    FMR1     X    146993447  146993679  232     0.0166    754           45332.5720
    ATXN10   22   46189527   46191972   2445    0.0000    0             0.0000

    
## fastxRepeatAnalysisReport.py
Generate "waterfall" and repeat count kde plots for unaligned CCS reads using flanking sequence.  
### Dependencies
 - [matplotlib](https://matplotlib.org/)
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [seaborn](https://seaborn.pydata.org/)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
### Usage
    $ python fastxRepeatAnalysisReport.py -h
    usage: fastxRepeatAnalysisReport.py [-h] [-t,--target TARGET]
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
                          FMR1,FUCHS,HTT,ALS,Sca10
      -m,--motifs MOTIFS  Comma separated motifs to search for. Will over-ride
                          presets. Default none.
      -o,--outDir OUTDIR  Output directory. Default cwd.
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
    m54006_180727_184845/10748847/ccs  2    27
    m54006_180727_184845/11469406/ccs  2    28
    m54006_180727_184845/11796868/ccs  2    27
    m54006_180727_184845/11928405/ccs  2    27
    m54006_180727_184845/12387293/ccs  1    326
    m54006_180727_184845/12649067/ccs  2    27
    m54006_180727_184845/12714233/ccs  2    27
    m54006_180727_184845/12780515/ccs  2    27
    m54006_180727_184845/13763208/ccs  2    27

[sample].[label].summary.csv

    $ column -ts, bc1002.FMR1.summary.csv
    totalReads     5679
    spanningReads  174
    oneSided       119
    poorAlignment  5386
    reference      resources/FMR1_L446_R503.fasta

## RepeatAnalysisReport.py
Generate "waterfall" and repeat count kde plots for aligned CCS reads using BED file with defined repeat regions.
### Dependencies
 - [matplotlib](https://matplotlib.org/)
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [seaborn](https://seaborn.pydata.org/)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
### Usage
    $ python RepeatAnalysisReport.py -h
    usage: RepeatAnalysisReport.py [-h] [-o,--outDir OUTDIR] [-s,--sample SAMPLE]
                                   [-f,--flanksize FLANKSIZE]
                                   inBAM inBED reference
    
    generate repeat-kde plot, waterfall plot, and summary for repeat expansions
    
    positional arguments:
      inBAM                 input BAM of CCS alignments
      inBED                 input BED defining location of repeat region. Repeats
                            ONLY, no flanking sequence. 
                            Columns (ordered): ctg, start, end, name, motifs
      reference             Reference fasta used for mapping BAM. Must have .fai
                            index.
    
    optional arguments:
      -h, --help            show this help message and exit
      -o,--outDir OUTDIR    Output directory. Default cwd.
      -s,--sample SAMPLE    Sample name, prepended to output files. Default None
      -f,--flanksize FLANKSIZE
                            Size of flanking sequence mapped for extracting repeat
                            region. Default 100
### Example BED file
    $ column -t resources/human_hs37d5.targets_repeatonly.bed
    4  3076554    3076717    HTT      CAG,CAA,CCG,CCA,CGG
    9  27573435   27573596   C9orf72  GGGGCC
    X  146993569  146993629  FMR1     CGG,AGG
### Output
_Same as fastxRepeatAnalysisReport.py_

## minimap2_e40.sh
Map repeat-extension data to targets, parameterized to reduce alignment cost for long expansions

### Dependencies
 - [minimap2](https://github.com/lh3/minimap2)
 - [samtools](https://github.com/samtools/samtools)

### Usage
    minimap2_e40.sh REFFASTA QUERY OUTBAM [SAMPLE]

## pbmm2_extension.sh
Map repeat-extension data to targets using pbmm2 (recommended), parameterized to reduce alignment cost for long expansions

### Dependencies
 - [pbmm2](https://github.com/PacificBiosciences/pbmm2) 

### Usage
    pbmm2_extension.sh REFFASTA QUERY OUTFILE [SAMPLE]

