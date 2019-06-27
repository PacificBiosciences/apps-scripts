# Description
This repo contains tools and wrapper scripts for processing, extracting and reporting sequence data generated with the PacBio No-Amp Targeted Sequencing Protocol (alpha).

Sequence read data generated with the No-Amp protocol on PacBio Sequel instruments make use of special asymmetric SMRTbell templates with different hairpin adapters on each end of the inserts and must be re-processed prior to analysis.  

This page describes the processing steps necessary and provides a set of extra tools for basic extraction and reporting of the results.  

Outputs from the analysis scripts include high-accuracy (>=QV20) CCS sequences for target regions so that users can easily analyze the results with other third party tools as necessary.   

## Dependencies

### Sequence Processing requires the following tools, available from [pbbioconda](https://github.com/PacificBiosciences/pbbioconda)
 - [recalladapters](https://github.com/pacificbiosciences/recalladapters/)
 - [lima](https://github.com/pacificbiosciences/barcoding)
 - [ccs](https://github.com/pacificbiosciences/unanimity/)
 - [pbmm2](https://github.com/PacificBiosciences/pbmm2/)

Additionally, the wrapper script described below makes use of **GNU parallel** to decrease time to results.

#### Notes
 - `recalladapters` 
 - `lima` and `ccs` are available in SMRT Link v6.0, but we recommend using the latest release, SMRT Link v7 or `pbccs` from the pbbioconda repository due to significant speed increases in recent versions
 - `pbmm2` is a minimap2 wrapper for native PacBio datasets (replacement for blasr)
 - The latest versions of lima,ccs, and pbmm2 tools will be available in SMRT Link v7

### Post-processing analysis tools are written in Python2.7 and require the following packages
 - [matplotlib](https://matplotlib.org/)
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [seaborn](https://seaborn.pydata.org/)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)

# Basic Repeat Analysis Workflow

## Raw Sequence Processing
Four steps are required to process sequence movies prior to repeat analysis:
1. Recall asymmetric SMRTbell adapters to properly separate subreads and allow for demultiplexing barcoded samples.
2. Demultiplex templates with barcodes on only one end of the insert.
   - This step can be skipped for samples without barcodes
3. Generate CCS reads and filter for reads >=QV20.
4. Align CCS reads to reference while allowing for gap extension in highly expanded alleles.

### Reference
For human samples, we recommend alignment to the reference `hs37d5`.  Example targets listed in the BED file below have `hs37d5` coordinates.  However, any reference will work so long as the BED file coordinates are paired with the reference used.

A wrapper shell script is provided to automate processing datasets prior to analysis.

    $ NoAmpRepeatPipeline.sh \
    movie.subreadset.xml \
    adapters.tetraloop.fasta \
    barcodes.fasta/.xml \
    reference.fasta/.mmi \
    output_directory

For datasets without barcodes, use the *no_Barcode* version

    $ NoAmpRepeatPipeline_noBarcode.sh \
    movie.subreadset.xml \
    adapters.tetraloop.fasta \
    reference.fasta/.mmi \
    output_directory

The NoAmp tetraloop adapter fasta can be found [here](resources/adapters.tetraloop.fasta), and the set of barcodes provided with the tetraloop adapters can be found [here](resources/barcodes_10.fasta) in the [resources subdirectory](resources/). 

### Output
After running the above wrapper scripts, results will be in four subdirectories under the provided output directory:
 - [output directory]
   - **align**   : contains the final aligned CCS reads for repeat analysis
   - **barcode** : contains unaligned, demuxed *subreads*, one per barcode used (only for barcoded data)
   - **ccs**     : contains unaligned CCS reads
   - **refarm**  : contains unaligned *subreads* and *scraps* with recalled adapters

## Repeat Analysis of CCS results
Once the data are processed as above, the aligned BAM files are used as inputs to the `RepeatAnalysisReport.py` script for extraction and reporting of results.  For barcoded data, a bash script `NoAmpRepeatReport.sh` is provided to loop over multiple samples.  

The script identifies reads spanning target regions and extracts the target region from each spanning read.  Motif counts per read are based on exact string matches.

    $ NoAmpRepeatReport.sh \
    output_directory \
    human_hs37d5.targets_repeatonly.bed \
    human_hs37d5.fasta

Where `output_directory` is the directory from the previous step. This script will add a subdirectory `reports` to the `output_directory`, and a further subdirectory for each barcoded sample

 - [output directory]
   - ... (from above)
   - reports
     - BC1007--BC1007
     - BC1016--BC1016
     - ...

Each barcoded directory contains results for one barcode and all targets listed in the BED file, including:
 - **summary.csv** : Table listing counts of reads spanning target regions
 - **repeatCounts.xlsx** : Excel file containing motif counts and total insert lengths, one sheet per target, one row per read.
 - **extractedSequence\_[target].fasta** : FASTA file of extracted inserts.  Name format: [movie]/[zmw]/ccs/[start]\_[stop]
 - __[target].*.png__ : Figures plotted from data in repeatCounts.xlsx 

### Un-barcoded samples
For data without barcodes, please use the python script directly (see below).

## Visualizing Repeats
We recommend [IGV v2.5.x](https://software.broadinstitute.org/software/igv/node/294) to visualize the alignments found in the `[output_directory]/align` folder.

## NEW SIMPLIFIED Repeat Analysis (beta)
Three new tools are provided to enable simplified reporting of repeat expansion.  
### Stream fastq records directly to waterfall script
    $ python extractRegion.py align/BC1026--BC1026.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python waterfall.py -m CGG,AGG -o FMR1.png

### Stream fastq directly to repeat counter
    $ python extractRegion.py align/BC1026--BC1026.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python countMotifs.py -m CGG,AGG | column -ts, | head
    readName                                   CGG  AGG  totalLength
    m54006_190116_193234/56885460/ccs/298_388  28   2    90
    m54006_190116_193234/17367745/ccs/299_391  29   1    92
    m54006_190116_193234/20840638/ccs/298_391  29   2    93
    m54006_190116_193234/9044396/ccs/300_393   29   2    93
    m54006_190116_193234/18481434/ccs/298_391  29   2    93
    m54006_190116_193234/18940755/ccs/298_391  29   2    93
    m54006_190116_193234/19136687/ccs/298_391  29   2    93
    m54006_190116_193234/20644190/ccs/298_391  29   2    93
    m54006_190116_193234/16778121/ccs/298_391  29   2    93

## extractRegion.py
This is a simple tool for extracting a target region from aligned CCS reads.  The target region is extracted by aligning the sequence immediately flanking the target to each CCS read intersecting the target.
### Dependencies
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
### Usage
    $ python extractRegion.py -h
    usage: extractRegion.py [-h] [-o,--outFq OUTFQ] [-f,--flanksize FLANKSIZE]
                            inBAM reference region
    
    extract target region from aligned BAMS using region flank alignments. Output
    format is fastq
    
    positional arguments:
      inBAM                 input BAM of CCS alignments
      reference             Reference fasta used for mapping BAM. Must have .fai
                            index.
      region                Target region, format '[chr]:[start]-[stop]'. Example
                            '4:3076604-3076660'
    
    optional arguments:
      -h, --help            show this help message and exit
      -o,--outFq OUTFQ      Output fastq file. Default stdout
      -f,--flanksize FLANKSIZE
                            Size of flanking sequence mapped for extracting repeat
                            region. Default 100

### Example
    $ python extractRegion.py m54006_190117_155211.refarm.barcoded.BC1026--BC1026.bam \
                              human_hs37d5.fasta \
                              4:3076604-3076660 | head
    @m54006_190117_155211/10616973/ccs/2259_2304
    CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG
    +
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @m54006_190117_155211/10879182/ccs/2259_2307
    CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG
    +
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @m54006_190117_155211/12386712/ccs/2259_2307
    CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG

## waterfall.py
This generates a waterfall-style plot from a fastx file (or stdin) of extracted repeat regions.

    $ python waterfall.py -h
    usage: waterfall.py [-h] [-i,--inFastx INFASTX] -o,--out OUT -m,--motifs
                        MOTIFS [-f,--format FORMAT] [-d,--dpi DPI]
    
    quick waterfall plot from fastx of extracted repeat sequences
    
    optional arguments:
      -h, --help            show this help message and exit
      -i,--inFastx INFASTX  Input Fastx file. Default stdin
      -o,--out OUT          Output file
      -m,--motifs MOTIFS    Search motifs, comma separated, most frequent first,
                            e.g. 'CGG,AGG'
      -f,--format FORMAT    Image format. Default png
      -d,--dpi DPI          Image resolution. Default 400

### Example
    $ python waterfall.py -i extractedSequence_FMR1.fasta -m CGG,AGG -o FMR1.png

## countMofifs.py
This plot will read a fastx file (or stdin) of extracted repeat regions and generate counts of exact motif string matches.

    $ python countMotifs.py -h
    usage: countMotifs.py [-h] [-i,--inFastx INFASTX] [-o,--out OUT] -m,--motifs
                          MOTIFS [-s,--sep SEP] [-r,--reverse]
    
    quick count of motifs per read from fastx of extracted repeat sequences
    
    optional arguments:
      -h, --help            show this help message and exit
      -i,--inFastx INFASTX  Input Fastx file. Default stdin
      -o,--out OUT          Output csv. Default stdout
      -m,--motifs MOTIFS    Search motifs, comma separated, most frequent first,
                            e.g. 'CGG,AGG'
      -s,--sep SEP          Field separator. Default ','
      -r,--reverse          Sort largest first. Default ascending order

### Example
    $ python countMotifs.py -i extractedSequence_FMR1.fasta -m CGG,AGG | column -ts,
    readName                                   CGG  AGG  totalLength
    m54006_190116_193234/56885460/ccs/298_388  28   2    90
    m54006_190116_193234/17367745/ccs/299_391  29   1    92
    m54006_190116_193234/20840638/ccs/298_391  29   2    93
    m54006_190116_193234/9044396/ccs/300_393   29   2    93
    m54006_190116_193234/18481434/ccs/298_391  29   2    93
    m54006_190116_193234/18940755/ccs/298_391  29   2    93
    m54006_190116_193234/19136687/ccs/298_391  29   2    93
    m54006_190116_193234/20644190/ccs/298_391  29   2    93
    m54006_190116_193234/16778121/ccs/298_391  29   2    93 

## RepeatAnalysisReport.py
Generate report scripts, extract repeat regions, plot "waterfall" and repeat count kde plots for aligned CCS reads using BED file with defined repeat region(s).
### Usage
    $ python RepeatAnalysisReport.py -h
    usage: RepeatAnalysisReport.py [-h] [-o,--outDir OUTDIR] [-s,--sample SAMPLE]
                                   [-f,--flanksize FLANKSIZE]
                                   inBAM inBED reference

    generate repeat-kde plot, waterfall plot, and summary for repeat expansions

    positional arguments:
      inBAM                 input BAM of CCS alignments
      inBED                 input BED defining location of repeat region. Repeats
                            ONLY, no flanking sequence. Columns: ctg, start, end,
                            name, motifs
      reference             Reference fasta used for mapping BAM. Must have .fai
                            index.

    optional arguments:
      -h, --help            show this help message and exit
      -o,--outDir OUTDIR    Output directory. Default cwd
      -s,--sample SAMPLE    Sample name, prepended to output files. Default None
      -f,--flanksize FLANKSIZE
                            Size of flanking sequence mapped for extracting repeat
                            region. Default 100

### Example BED file
    $ column -t resources/human_hs37d5.targets_repeatonly.bed
    4   3076604    3076660    HTT      CAG,CAA,CCG,CCA,CGG
    9   27573435   27573596   C9orf72  GGGGCC
    X   146993569  146993628  FMR1     CGG,AGG
    22  46191235   46191304   ATXN10   ATTCT,ATTCC,ATTTCT,ATTCCT

Note that the **first** motif listed for each target is the primary motif.  All subsequent motifs are potential interruption motifs within the repeat region.  Start/Stop locations are 1-based (as in IGV) and **_do not_** include flanking sequence.  

__*Columns:*__ **contig**, **start**, **stop**, **name**, **motifs**

### Output

[SampleBarcode].extractedSequence\_[Target].fasta
 - One FASTA file for each target in the BED file.  This is the part of the sequence read that corresponds to the target region for each read.
 - Naming convention: [movie]/[zmw]/ccs/[start]\_[stop]

[SampleBarcode].[Target].insertSize_dist.png

![Target Region Size](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/examples/outputs/ExampleBC.FMR1.insertSize_dist.png)

[SampleBarcode].[Target].motifCount_dist.png

![Target Region Motif Counts](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/examples/outputs/ExampleBC.FMR1.motifCount_dist.png)

[SampleBarcode].[Target].waterfall.png

![Waterfall plot Example](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/examples/outputs/ExampleBC.FMR1.waterfall.png)

[sample].[label].repeatCounts.csv

![Waterfall plot Example](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/examples/outputs/ExampleBC.repeatCounts.xlsx_screenshot.png)

[SampleBarcode].summary.csv

    $ column -ts, ExampleBC.summary.csv  
    target   totalReads  oneSided  poorAlignment  spanningReads
    ATXN10   112         0         0              112
    C9orf72  0           0         0              0
    FMR1     72          0         0              72
    HTT      187         0         0              187

# Additional Stand-alone Tools
The following are provided as convenience tools



## NoAmpRestrictionDiagnostics.py
This tool can be used to help diagnose issues in the background-reducing digest.  See example [restriction_enzymes.tsv](resources/restriction_enzymes.tsv) in `resources` folder.  Note that enzymes with degenerate cutsites are listed multiple times.
### Dependencies
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
 - [pbcore](https://github.com/pacificbiosciences/pbcore/)
### Usage
    $ python NoAmpRestrictionDiagnostics.py  -h
    usage: NoAmpRestrictionDiagnostics.py [-h] [-o,--outDir OUTDIR]
                                          [-s,--subreadsBAM SUBREADSBAM]
                                          [-a,--adapterFasta ADAPTERFASTA]
                                          [-m,--minAdapterScore MINADAPTERSCORE]
                                          [-f,--minFlankScore MINFLANKSCORE]
                                          inBAM inBED reference reTable
    
    Generate restriction enzyme and [optionally] adapter report for NoAmp protocol
    
    positional arguments:
      inBAM                 input BAM of CCS alignments to reference
      inBED                 input BED defining location of repeat region(s).
                            Columns ['ctg','start','end','name']
      reference             Reference fasta used for mapping BAM
      reTable               tsv table of restriction enzymes and cut sequences.
                            Use multiple lines for degenerate cutsites Columns:
                            ['name','cutsite']
    
    optional arguments:
      -h, --help            show this help message and exit
      -o,--outDir OUTDIR    Output directory. Default cwd
      -s,--subreadsBAM SUBREADSBAM
                            refarmed subreads bam for identifying adapter types.
                            Must have 'ad' tag from '--adpqc' option in bam2bam.
                            Default None
      -a,--adapterFasta ADAPTERFASTA
                            fasta of re-farmed adapter sequences. Required if '-s'
                            option used. Default None
      -m,--minAdapterScore MINADAPTERSCORE
                            minimum adapter score. Default 0
      -f,--minFlankScore MINFLANKSCORE
                            minimum adapter flanking score. Default 0
### Example
    $ python NoAmpRestrictionDiagnostics.py ccs.recalled.2.hs37d5.bam \
                                            human_hs37d5.targets.bed \
                                            human_hs37d5.fasta \
                                            restriction_enzymes.tsv \
                                            -s tetraloop_recall.subreads.bam \
                                            -a adapters.tetraloop.fasta
### Output
    $ column -t adapterReport.tsv
    adapters     ZMWcount  TotalFraction  mappedCCS  CCSfrac
    NoAd;NoAd    96376     74.96%                    0.00%
    TL_42;tc6    16556     12.88%         11277.00   68.11%
    NoAd;TL_42   11097     8.63%          1.00       0.01%
    NoAd;tc6     3122      2.43%                     0.00%
    tc6;tc6      816       0.63%          441.00     54.04%
    TL_42;TL_42  607       0.47%          361.00     59.47%

    $ column -t restrictionCounts.csv
    target  adapters     Acc1_int  SexA1_int  BamH1_int  EcoR1_int  EcoR1_rend  ZMWcount  Fraction
    ATXN10  TL_42;tc6    2         1          1          0          0           1         0.06%
    ATXN10  TL_42;tc6    2         1          1          0          1           434       24.00%
    ATXN10  tc6;tc6      2         1          0          0          0           1         0.06%
    ATXN10  tc6;tc6      2         1          1          1          0           3         0.17%
    FMR1    TL_42;TL_42  0         0          1          0          1           1         0.06%
    FMR1    TL_42;tc6    0         0          1          0          0           1         0.06%
    FMR1    TL_42;tc6    0         0          1          0          1           189       10.45%
    FMR1    TL_42;tc6    0         0          1          1          0           140       7.74%
    FMR1    tc6;tc6      0         0          1          1          0           3         0.17%
    HTT     TL_42;tc6    0         0          1          0          0           10        0.55%
    HTT     TL_42;tc6    0         0          1          0          1           602       33.30%
    HTT     TL_42;tc6    0         0          1          1          0           421       23.29%
    HTT     tc6;tc6      0         0          1          0          1           2         0.11%

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
    HTT      4    3076554    3076717    163     0.0167    1164          69787.1946
    C9orf72  9    27573435   27573596   161     0.0167    0             0.0000
    FMR1     X    146993447  146993679  232     0.0161    754           46761.0800
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

## pbmm2_extension.sh
Map repeat-extension data to targets using pbmm2 (recommended), parameterized to reduce alignment cost for long expansions

### Dependencies
 - [pbmm2](https://github.com/PacificBiosciences/pbmm2) 

### Usage
    pbmm2_extension.sh REFFASTA QUERY OUTFILE [SAMPLE]

## minimap2_e40.sh
Map repeat-extension data to targets, parameterized to reduce alignment cost for long expansions

### Dependencies
 - [minimap2](https://github.com/lh3/minimap2)
 - [samtools](https://github.com/samtools/samtools)

### Usage
    minimap2_e40.sh REFFASTA QUERY OUTBAM [SAMPLE]

