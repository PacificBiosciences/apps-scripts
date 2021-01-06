# Description
This repository contains instructions for processing and repeat analysis of sequence data generated with the PacBio No-Amp Targeted Sequencing Protocol with simplified double Cas9 cut.

UPDATE: RepeatAnalysis Tools in this repository now use Python 3. 

Outputs from the analysis scripts include high-accuracy (>=QV20) CCS sequences for target regions so that users can easily analyze the results with other third party tools as necessary.   

## Dependencies

### Sequence Processing on the CL requires the following tools, available from [pbbioconda](https://github.com/PacificBiosciences/pbbioconda)
 - [ccs](https://ccs.how/)
 - [lima](https://github.com/pacificbiosciences/barcoding)
 - [pbmm2](https://github.com/PacificBiosciences/pbmm2/)

### UPDATE Repeat analysis tools are written in Python 3 and require the following packages
 - [matplotlib](https://matplotlib.org/)
 - [numpy](http://www.numpy.org/)
 - [pandas](https://pandas.pydata.org/)
 - [seaborn](https://seaborn.pydata.org/)
 - [mappy](https://github.com/lh3/minimap2/tree/master/python)
 - [pysam](https://pysam.readthedocs.io/en/latest/index.html)
 - [scikit-learn](https://scikit-learn.org/stable/index.html)
 - [pbcore](https://github.com/pacificbiosciences/pbcore/)

# Example Dataset
An example dataset with repeat-expansion genotypes for HTT and FMR1 can be found [here](https://downloads.pacbcloud.com/public/dataset/RepeatExpansionDisorders_NoAmp/)

# Basic Repeat Analysis Workflow

## Raw Sequence Processing
Repeat analysis requires 3 basic steps
1. CCS (HiFi)
2. Demultiplexing
3. Alignment

We recommend using the latest [CCS](https://ccs.how/) version for processing repeat-expansion sequence on the command line. 

With improvements to CCS versions 5+, previously problematic very long repeat expansions now generate output ccs reads on default settings.  With the addition of the `--all` option to CCS, any sequencing well that produces subreads will have representative data in the output of ccs.  See [CCS --all](https://ccs.how/faq/mode-all.html) for details. We are continually trying to improve on the CCS algorithm, so please contact us if you encounter cases where you see many unpolished consensus reads in ccs outputs. 

For customers still using CCS v4.2 (part of SMRT Link v9.0.0), we recommend upgrading from the links above.  If you must use v4.2, please add the options noted below.

The wrapper script `preprocess.sh` is provided below as a single-command for running all three preprocessing steps.

### CCS
In order to get the maximum yield for the longest expansion alleles, we need to run ccs with heuristics turned off and full-length draft mode:

    $ ccs in.subreads.bam out.ccs.bam --all

For version 4.2

    $ ccs in.subreads.bam out.ccs.bam --disable-heuristics --draft-mode full

Running ccs in this mode will be considerably slower than running on default, so it is recommended to chunk the ccs calls if you have access to a cluster. The ccs program has built-in chunking capabilities. See the script `chunkCCS.sh` for an example of chunking ccs on the command line.

This command will start 16 chunked jobs with 9 threads each using qsub: 

    $ ./chunkCCS.sh 16 9 cluster ccs in.subreads.bam out.ccs.bam --all

For version 4.2

    $ ./chunkCCS.sh 16 9 cluster ccs in.subreads.bam out.ccs.bam --disable-heuristics --draft-mode full

This command will start 4 ccs jobs with 4 threads each on the local machine

    $ ./chunkCCS.sh 4 4 local ccs in.subreads.bam out.ccs.bam --all

For version 4.2

    $ ./chunkCCS.sh 4 4 local ccs in.subreads.bam out.ccs.bam --disable-heuristics --draft-mode full

### Demultiplex
Demultiplexing of the CCS reads is done with the program `lima`.

    $ lima --ccs --same                   \
           --split-bam-named --peek-guess  \
           inCCS.bam barcodes.fasta demuxCCS.bam

A summary of the demultiplexing outputs can be found in the `lima` output directory with the suffix \*summary. 

### Align to Reference
In order to align long expansion repeats through short reference seqeunce, we provide a set of modified parameters to the program `pbmm2` for alignment.

    $ pbmm2_extention.sh human_hs37d5.fasta ccs.bam mapped.ccs.bam sampleName

### Reference
For human samples, we recommend alignment to the reference `hs37d5`.  Example targets listed in the BED file below have `hs37d5` coordinates.  However, any reference will work so long as the target coordinates are paired with the reference used.

### Preprocess Script
A wrapper script `preprocess.sh` is provided for one-step preprocessing from subreads -> ccs -> demux -> align.  This script can run "local" on the executing machine or will submit to a cluster using qsub.

To run preproccesing on the cluster with 32 chunks and 8 procs per chunk:

    $ preprocess.sh subreads.bam|xml barcodes.fasta|xml reference.fasta outdir 32 8 cluster

Following completion of the script, three folders will be created in the output directory `ccs`, `demux`, and `align`.  BAM files in the `align` subdirectory are ready for analysis in the repeat analysis reporting scripts below.

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
    HTT      4    3076604    3076693    89      0.0018    408           227454.0622 
    C9orf72  9    27573435   27573596   161     0.0017    210           120320.2169 
    FMR1     X    146993569  146993628  59      0.0018    242           135015.6210 
    ATXN10   22   46191235   46191304   69      0.0018    268           149907.1604 
    
### Example BED file
    $ column -t resources/human_hs37d5.targets_repeatonly.bed
    4   3076604    3076693    HTT      CAG,CAA,CCG,CCA,CGG          0
    9   27573435   27573596   C9orf72  GGGGCC,GGGGCG                1
    X   146993569  146993628  FMR1     CGG,AGG                      0
    22  46191235   46191304   ATXN10   ATTCT,ATTCCT,ATTTCT,ATTCC    0

Columns are: {chr} {start} {stop} {name} {motifs} {revcomp}

The final column `revcomp` is optional and indicates whether the aligned sequence should be reverse complemented (1) or not (0) before exporting extracted sequence.  This is useful in cases where the reference sequence used is the opposite strand compared to the typical orientation.  

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

## clusterByRegion.py
K-means clustering of reads based on kmer counts over the repeat region of interest provides a reliable way to phase alleles.  Output includes haplotagged BAM (tag="HP") and summary stats for target motifs.
### Usage
    $ python clusterByRegion.py -h
    usage: clusterByRegion.py [-h] -m,--motifs MOTIFS [-k,--kmer KMER]
                          [-c,--clusters CLUSTERS] [-r,--revcomp]
                          [-p,--prefix PREFIX] [-f,--flanksize FLANKSIZE]
                          [-s,--seed SEED] [-x,--noBam] [-d,--drop]
                          [-u,--collapseHP]
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
      -r,--revcomp          Reverse complement extracted sequence
      -p,--prefix PREFIX    Output prefix. Default ./cluster
      -f,--flanksize FLANKSIZE
                            Size of flanking sequence mapped for extracting repeat
                            region. Default 100
      -s,--seed SEED        Seed for resampling ci95. Default 42
      -x,--noBam            Do not export HP-tagged bam of clustered reads
      -d,--drop             Drop reads with no cluster in output bam. Default keep
                            all reads.
      -u,--collapseHP       Collapse homopolymers before analysis. Default use
                            original sequence.

### Example
    $ python clusterByRegion.py -m CGG,AGG \
                                -p cluster/FMR1 \
                                -d \
                                combined.consensusalignmentset.bam \
                                human_hs37d5.fasta \
                                'X:146993569-146993628'

    $ column -ts, cluster/FMR1.summary.csv 
                          CGG     CGG    CGG          AGG     AGG   AGG      totalBp  totalBp  totalBp       Read
                          median  mean   ci95         median  mean  ci95     median   mean     ci95          count
    cluster0_numreads86   320.5   320.1  (266 - 384)  1       1.1   (0 - 3)  985.5    988.0    (857 - 1168)  86
    cluster1_numreads151  27.0    27.0   (26 - 28)    2       2.0   (2 - 2)  87.0     87.2     (84 - 90)     151
    
    $ head cluster/FMR1.readnames.txt
    >cluster0_numreads86
    m64012_190806_011308/22087479/ccs/295_1272
    m64012_190806_011308/172818496/ccs/294_1254
    m64012_190806_011308/48760706/ccs/296_1331
    m64012_190806_011308/180027489/ccs/289_1286
    m64012_190806_011308/43582972/ccs/295_1228


![Haplotagged Bam](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/bc1019.FMR1_haptagged.png)

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
                            [-r,--revcomp]
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
      -r,--revcomp          Rev-comp extracted region. Default Reference Direction

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
