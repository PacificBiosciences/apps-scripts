# Scripts for manipulating FALCON and FALCON-Unzip Assemblies

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
    
## Align haplotigs to their primary contig with nucmer

### Dependencies
[samtools](http://samtools.sourceforge.net/)

[mummer](http://mummer.sourceforge.net/)

### Usage
    alignHaplotigs2Primary.sh myUnzipAsm.fasta primaryContigID

primary contig ID formatted as: 000123F

### Visualization

view filtered delta file (primaryContigID.g.delta) in [assemblytics](http://qb.cshl.edu/assemblytics/) or [mummerplot](http://mummer.sourceforge.net/manual/#mummerplot).
