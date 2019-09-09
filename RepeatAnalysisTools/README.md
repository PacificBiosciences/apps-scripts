# Description
This repository contains instructions for processing and repeat analysis of sequence data generated with the PacBio No-Amp Targeted Sequencing Protocol with simplified double Cas9 cut.

Outputs from the analysis scripts include high-accuracy (>=QV20) CCS sequences for target regions so that users can easily analyze the results with other third party tools as necessary.   

## Dependencies

SMRT-Link v7.0+ is required for alignment of expanded allele sequences.  This and later versions make use of pbmm2 (PacBio minimap2 wrapper) for alignment.  Viewing of aligned BAMS is improved via parameter over-rides to accomodate mapping long extended alleles through short reference targets.   

### Sequence Processing can also be accomplished on the CL using the following tools, available from [pbbioconda](https://github.com/PacificBiosciences/pbbioconda)
 - [lima](https://github.com/pacificbiosciences/barcoding)
 - [ccs](https://github.com/pacificbiosciences/unanimity/)
 - [pbmm2](https://github.com/PacificBiosciences/pbmm2/)

### Repeat analysis tools are written in Python2.7 and require the following packages
 - [matplotlib](https://matplotlib.org/)
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [seaborn](https://seaborn.pydata.org/)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
 - [scikit-learn](https://scikit-learn.org/stable/index.html)

# Basic Repeat Analysis Workflow

## Raw Sequence Processing
Repeat analysis requires 3 basic steps in SMRT-Link
1. Demultiplex
2. CCS
3. Alignment

### Reference
For human samples, we recommend alignment to the reference `hs37d5`.  Example targets listed in the BED file below have `hs37d5` coordinates.  However, any reference will work so long as the target coordinates are paired with the reference used.

## Demultiplexing in SMRT-Link
1. Select 'Create New Analysis'
2. Select your raw data and give the job a name (Next)
3. Select 'Demultiplex Barcodes' and the proper barcode set
4. Accept default parameters (Same on both sides, Infer barcodes)

![Demultiplex Barcodes](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/demux_jobSetup.png)

Check the results once the job is completed.

![Demultiplex Results](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/demuxResults.png)

## CCS and Mapping
1. Select 'Create New Analysis'
2. Select demuxed data and give the job a name.  
    * *Note* Select each barcoded sample from the child menu after clicking the number in the field 'Demultiplexed Subsets'.
![Select Each Demuxed Set](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/ccsDataSelect2.png)
3. Select the option 'One Analysis per Data Set -- Identical Parameters' (Next)
4. Select 'CCS with Mapping' and the desired reference for alignment
5. Select 'Advanced Analysis Parameters'
    * Turn 'Consolidate BAM' to *ON* to generate a single aligned BAM for output.
![Consolidate BAM](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/ccsParams1.png)
    * Set Minimum Concordance to 0
    * Add override options '-L 0.1 -E 0'
![Alignment Parameter Change](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/ccsParams2.png)


## Repeat Analysis of CCS results
Once the data are processed as above, the aligned BAM files are used as inputs to the following scripts for extraction and reporting of results.

### Visualizing Repeats
We recommend [IGV v2.5.x](https://software.broadinstitute.org/software/igv/node/294) (or later) to visualize the aligned BAM file generated for each CCS Mapping job.

![Visualization in IGV](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/c9orf72_IGV_med.png)
Repeat expansions are clearly labeled in IGV, and reads can be grouped into haplotypes.

## Auto-Generate All Reports
The script `makeReports.sh` is provided to generate all reports in a single command.  Please see below for examples of individual tools as well as target BED format.  This script uses the default `python` instance in the command path -- edit the script to use a different instance as needed.

### Usage
    $ makeReports.sh -h

    Usage: makeReports.sh targets.BED reference.fasta outputDir aligned1.bam [aligned2.bam aligned3.bam ...] 

The file `[outputDir]/runReportCmds.sh` contains commands used for generating the reports.

## clusterByRegion.py
NEW! K-means clustering of reads based on kmer counts over the repeat region of interest provides a reliable way to phase alleles.  Output includes haplotagged BAM (tag="HP") and summary stats for target motifs.
### Usage
    $ python clusterByRegion.py -h
    usage: clusterByRegion.py [-h] -m,--motifs MOTIFS [-k,--kmer KMER]
                              [-c,--clusters CLUSTERS] [-p, --prefix PREFIX]
                              [-f,--flanksize FLANKSIZE] [-x,--noBam] [-d,--drop]
                              inBAM reference region
    
    kmer clustering by target region
    
    positional arguments:
      inBAM                 input BAM of CCS alignments
      reference             Reference fasta used for mapping BAM. Must have .fai
                            index.
      region                Target region, format '[chr]:[start]-[stop]'. Example
                            '4:3076604-3076660'
    
    optional arguments:
      -h, --help            show this help message and exit
      -m,--motifs MOTIFS    comma-separated list of motifs to count
      -k,--kmer KMER        kmer size for clustering. Default 3
      -c,--clusters CLUSTERS
                            clusters/ploidy count. Default 2
      -p, --prefix PREFIX   Output fastq file. Default ./cluster
      -f,--flanksize FLANKSIZE
                            Size of flanking sequence mapped for extracting repeat
                            region. Default 100
      -x,--noBam            Do not export HP-tagged bam of clustered reads
      -d,--drop             Drop reads with no cluster in output bam. Default keep
                            all reads.

### Example
    $ python clusterByRegion.py -m CGG,AGG \
                                -p cluster/FMR1 \
                                -d \
                                combined.consensusalignmentset.bam \
                                human_hs37d5.fasta \
                                'X:146993569-146993628'
    $ column -ts, cluster/bc1019--bc1019.mapped.FMR1.summary.csv 
                          CGG     CGG    CGG          AGG     AGG   AGG      totalBp  totalBp  totalBp       Read
                          median  mean   ci95         median  mean  ci95     median   mean     ci95          count
    cluster0_numreads86   320.5   320.1  (266 - 384)  1       1.1   (0 - 3)  985.5    988.0    (857 - 1168)  86
    cluster1_numreads151  27.0    27.0   (26 - 28)    2       2.0   (2 - 2)  87.0     87.2     (84 - 90)     151


![Haplotagged Bam](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/bc1019.FMR1_haptagged.png)


## countOnTarget.py
Generate table of ZMW counts per target.
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
    HTT      4    3076604    3076660    56      0.0018    408           227454.0622 
    C9orf72  9    27573435   27573596   161     0.0017    210           120320.2169 
    FMR1     X    146993569  146993628  59      0.0018    242           135015.6210 
    ATXN10   22   46191235   46191304   69      0.0018    268           149907.1604 
    
### Example BED file
    $ column -t resources/human_hs37d5.targets_repeatonly.bed
    4   3076604    3076660    HTT      CAG,CAA,CCG,CCA,CGG
    9   27573435   27573596   C9orf72  GGGGCC
    X   146993569  146993628  FMR1     CGG,AGG
    22  46191235   46191304   ATXN10   ATTCT,ATTCC,ATTTCT,ATTCCT

## coveragePlot.py
Generate coverage plot of results. (Example has 4 multiplexed targets)
### Example
    $python coveragePlot.py combined.consensusalignmentset.bam 

![Coverage Plot](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/4target.coveragePlot.png)

Add target labels:

    $python coveragePlot.py combined.consensusalignmentset.bam -t resources/human_hs37d5.targets.bed

![Coverage Plot](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/4target.coveragePlot_wlabels.png)

##
The following tools are provided to enable simplified reporting of repeat expansion.  
Click here for [Previous Repeat Reporting Tools](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/previous/).

### Extract fastq of repeat regions.
    $ python extractRegion.py combined.consensusalignmentset.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' | head
    @m64012_190806_011308/66977880/ccs/296_380
    CGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGG
    +
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~~~~~~~~~~z~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @m64012_190806_011308/178717562/ccs/296_383
    CGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGG
    +
    ~~~~~~~~~y~~~~~~~~u~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @m64012_190806_011308/39650860/ccs/295_382
    CGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGAGGCGGCGGCGGCGGCGGCGGCGGCGGCGG

### Stream fastq records directly to waterfall script
    $ python extractRegion.py combined.consensusalignmentset.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python waterfall.py -m CGG,AGG -o FMR1.png

![WaterFall Plot](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/FMR1.waterfall.png)

### Stream fastq directly to repeat counter
    $ python extractRegion.py combined.consensusalignmentset.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python countMotifs.py -m CGG,AGG | column -ts, 
    readName                                     CGG  AGG  totalLength
    m64012_190806_011308/39650860/ccs/295_382    27   2    87
    m64012_190806_011308/39651811/ccs/295_382    27   2    87
    m64012_190806_011308/41092035/ccs/295_382    27   2    87
    m64012_190806_011308/155584458/ccs/296_383   27   2    87
    m64012_190806_011308/46597674/ccs/295_382    27   2    87
    ...
    m64012_190806_011308/169740995/ccs/294_1298  326  1    1004
    m64012_190806_011308/131729381/ccs/295_1300  333  1    1005
    m64012_190806_011308/99680942/ccs/294_1299   323  1    1005
    m64012_190806_011308/134742966/ccs/296_1302  331  1    1006
    m64012_190806_011308/177735643/ccs/282_1290  328  1    1008
    m64012_190806_011308/8391101/ccs/295_1304    333  1    1009

### Stream fastq directly to count plotter
    $ python extractRegion.py combined.consensusalignmentset.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python plotCounts.py -m CGG -n FMR1

![hist.motifCount.png](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/hist.motifCount.png)
![hist.insertSize.png](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/hist.insertSize.png)

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
    $ python extractRegion.py combined.consensusalignmentset.bam \
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
### Dependencies
 - [pbcore](https://github.com/PacificBiosciences/pbcore)
 - [matplotlib](https://matplotlib.org/3.1.0/users/installing.html)
 - [numpy](https://www.numpy.org/)

### Usage
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

## countMotifs.py
This plot will read a fastx file (or stdin) of extracted repeat regions and generate counts of exact motif string matches.
### Dependencies
 - [pbcore](https://github.com/PacificBiosciences/pbcore)
 - [numpy](https://www.numpy.org/)

### Usage
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

## plotCounts.py
Plot histograms of repeated motif counts and total repeat region
### Dependencies
 - [pbcore](https://github.com/PacificBiosciences/pbcore)
 - [matplotlib](https://matplotlib.org/3.1.0/users/installing.html)
 - [seaborn](https://seaborn.pydata.org/)

### Usage
    $ python plotCounts.py -h
    usage: plotCounts.py [-h] [-i,--inFastx INFASTX] [-o,--out OUT] -m,--motif
                         MOTIF [-n,--name NAME] [-f,--format FORMAT]
                         [-d,--dpi DPI]
    
    generate histograms of motif counts and expansion size
    
    optional arguments:
      -h, --help            show this help message and exit
      -i,--inFastx INFASTX  Input Fastx file. Default stdin
      -o,--out OUT          Output prefix. default 'hist'
      -m,--motif MOTIF      Search motif
      -n,--name NAME        Title/name for figure
      -f,--format FORMAT    Image format. Default png
      -d,--dpi DPI          Image resolution. Default 400

### Example
    $ python plotCounts.py -m CGG -n FMR1 -i extractedSequence_FMR1.fasta
