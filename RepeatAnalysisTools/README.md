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

### Un-barcoded samples
For data without barcodes, please use the python script directly (see below).

## Repeat Analysis of CCS results
Once the data are processed as above, the aligned BAM files are used as inputs to the following scripts for extraction and reporting of results.

### Visualizing Repeats
We recommend [IGV v2.5.x](https://software.broadinstitute.org/software/igv/node/294) to visualize the alignments found in the `[output_directory]/align` folder.


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
    
### Example BED file
    $ column -t resources/human_hs37d5.targets_repeatonly.bed
    4   3076604    3076660    HTT      CAG,CAA,CCG,CCA,CGG
    9   27573435   27573596   C9orf72  GGGGCC
    X   146993569  146993628  FMR1     CGG,AGG
    22  46191235   46191304   ATXN10   ATTCT,ATTCC,ATTTCT,ATTCCT

Note that the **first** motif listed for each target is the primary motif.  All subsequent motifs are potential interruption motifs within the repeat region.  Start/Stop locations are 1-based (as in IGV) and **_do not_** include flanking sequence.  

__*Columns:*__ **contig**, **start**, **stop**, **name**, **motifs**

## NEW SIMPLIFIED Repeat Analysis (beta)
Three new tools are provided to enable simplified reporting of repeat expansion.  
Click here for [Previous Repeat Reporting Tool](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/previous/).

### Extract fastq of repeat regions.
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

### Stream fastq records directly to waterfall script
    $ python extractRegion.py align/BC1026--BC1026.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python waterfall.py -m CGG,AGG -o FMR1.png

![WaterFall Plot](https://github.com/PacificBiosciences/apps-scripts/blob/master/RepeatAnalysisTools/images/FMR1.waterfall.png)

### Stream fastq directly to repeat counter
    $ python extractRegion.py align/BC1026--BC1026.bam \
                              human_hs37d5.fasta \
                              'X:146993569-146993628' \
                              | python countMotifs.py -m CGG,AGG | column -ts, 
    readName                                   CGG  AGG  totalLength
    m54006_190116_193234/56885460/ccs/298_388  28   2    90
    m54006_190116_193234/17367745/ccs/299_391  29   1    92
    m54006_190116_193234/20840638/ccs/298_391  29   2    93
    m54006_190116_193234/9044396/ccs/300_393   29   2    93
    m54006_190116_193234/18481434/ccs/298_391  29   2    93
    m54006_190116_193234/16778121/ccs/298_391  29   2    93
    ...
    m54006_190116_193234/21430836/ccs/296_464  56   0    168
    m54006_190116_193234/19661344/ccs/295_463  56   0    168
    m54006_190116_193234/14418501/ccs/295_463  56   0    168
    m54006_190116_193234/42009135/ccs/295_463  56   0    168
    m54006_190116_193234/62587597/ccs/294_462  56   0    168
    m54006_190116_193234/56754595/ccs/295_463  56   0    168

### Stream fastq directly to count plotter
    $ python extractRegion.py align/BC1026--BC1026.bam \
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
### Dependencies
 - [pbcore](https://github.com/PacificBiosciences/pbcore)
 - [matplotlib](https://matplotlib.org/3.1.0/users/installing.html)
 - [numpy](https://www.numpy.org/)

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

### Example
    $ python plotCounts.py -m CGG -n FMR1 -i extractedSequence_FMR1.fasta
