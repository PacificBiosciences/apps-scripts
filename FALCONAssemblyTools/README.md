# SCRIPTS FOR MANIPULATING FALCON AND FALCON-UNZIP ASSEMBLIES


## Split Primary Contigs and Haplotigs into Separate Fasta Files

### Dependencies
[samtools](http://samtools.sourceforge.net/)

### Usage

    splitPrimariesHaplotigs.sh myUnzipAsm.fasta




## Filter Contigs with Poor Polishing

### Dependencies
[pbcore](https://github.com/PacificBiosciences/pbcore)

Default setting is to remove contigs with fewer than 50% of their bases polished. During arrow polishing, base positions with fewer than 5 raw reads mapped will not have an updated consensus base called, but rather, will have the reference base returned in lowercase. This script filters out contigs that have low raw read support.

### Usage

    trimLowercaseContigs.py myPolishedAsm.fasta myoutfile.fasta
    
    
    
    
## Align Haplotigs/Associated Contigs to their Primary Contig with Nucmer

### Dependencies
[samtools](http://samtools.sourceforge.net/)

[mummer](http://mummer.sourceforge.net/)

### Usage
For Unzip Assemblies

    alignHaplotigs2Primary.sh myUnzipAsm.fasta primaryContigID

For FALCON Assemblies

    alignAssoc2Primary.sh myFALCONAsm.fasta primaryContigID


Assembly file must contain both primary contigs and haplotigs or associated contigs.

Primary contig ID formatted as: 000123F.

### Visualization

view filtered delta file (primaryContigID.g.delta) in [assemblytics](http://qb.cshl.edu/assemblytics/) or [mummerplot](http://mummer.sourceforge.net/manual/#mummerplot).




## Remove Nested Haplotigs

### Dependencies
python2

[numpy](http://www.numpy.org/)

[samtools](http://samtools.sourceforge.net/)

[mummer](http://mummer.sourceforge.net/)

### Usage
    removeNestedHaplotigs.sh myUnzipAsm.fasta

This shell script aligns all haplotigs against their primary contig using nucmer, then identifies haplotigs that are nested within other haplotigs using nestedHaplotigs.py. Understanding why these small haplotigs are assembled is an area of ongoing research, but users may wish to remove these haplotigs from their assembly.




## Mean Contig Qualities
Outputs contig ID, mean quality across all bases, and contig length

### Usage 
    fastq2meanQ.py myAsm.fastq

Fastq sequence and quality lines must not be wrapped



## Identify Homologous contigs in primary contigs output

### Depedencies
python2.7.x

[mummer 4.0.0](https://github.com/mummer4/mummer/releases/tag/v4.0.0beta)

### Usage
    python get_homologs.py --nproc 24 /path/to/primary_contigs.fasta

The output is a directory of mummerplot 'Qfiles' listing homologous contigs, along with a plots.sh script that
can be run independently to draw dotplots for each reference contig with all homologous queries aligned to it.