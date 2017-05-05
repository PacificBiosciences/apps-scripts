# A collection of scripts and tools for manupulating FALCON and FALCON-Unzip Assemblies

## Split Primary Contigs and Haplotigs into Separate Fasta Files

## Dependencies
[samtools](http://samtools.sourceforge.net/)


    splitPrimariesHaplotigs.sh myUnzipAsm.fasta

## Filter contigs with poor polishing

## Dependencies
[pbcore](https://github.com/PacificBiosciences/pbcore)

Default setting is to remove contigs with fewer than 50% of their bases polished. During arrow polishing, base positions with fewer than 5 raw reads mapped will not have an updated consensus base called, but rather, will have the reference base returned in lowercase. This script filters out contigs that have low raw read support.

    trimLowercaseContigs.py myPolishedAsm.fasta myPolishedFilteredAsm.fasta
