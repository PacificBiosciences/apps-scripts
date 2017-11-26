# Reseq_BigGenome
## Resequencing hack for draft genomes larger than 4.29 Gb

The resequencing pipeline is recommended for polishing genomes assembled with PacBio data to increase base accuracy. Resequencing uses [BLASR](https://github.com/PacificBiosciences/blasr) to map raw reads to the draft reference and arrow for [genomic consensus](https://github.com/PacificBiosciences/GenomicConsensus) and can be accessed with the [SMRT Link GUI](http://www.pacb.com/wp-content/uploads/SMRT_Link_User_Guide.pdf) or [pbsmrtpipe](http://pbsmrtpipe.readthedocs.io/en/master/getting_started.html#basic-resequencing) at the command line. BLASR cannot index reference sequences longer than 2^32 (4.29 Gb).

This repository contains an unsupported workflow which will help you split your reference into two set of contigs, run resequencing on the two contig sets, and then merge the resequencing jobs back together by removing duplicate alignments. Specficially, reads that map to both contig sets are identified, and only the better alignment is retained. From this stage, you can run genomic consensus, or estimate coverage across your assembly, among other things. I have made the final merged BAM file optional but have included some code which can be uncommented.

### A note on resequencing
If you are using this pipeline, you have a large genome and a lot of raw reads. You can imagine that finding reads that map to both contig sets ("shared reads") can be inefficient. Fortunately, the engineering of read alignment in the resequencing pipeline allow the process of finding and assessing the alignments of shared reads to be done in a parallel manner.

During raw read alignment, the resequencing pipeline chunks subreads into smaller datasets and runs pbalign (BLASR) in parallel on each dataset chunk. This stage produces a series of files: `<job_dir>/tasks/pbalign.tasks.pbalign-[0-9]{1,2}/mapped.alignmentset.bam`. The subread chunking process is reproducible such that, for example, the file `pbalign.tasks.pbalign-13/mapped.alignmentset.bam` from resequencing of contig set1 and contig set 2 contain the same subreads. Thus, the process of finding shared read alignments and ommiting the poorer alignment can be done in parallel on each pair of `mapped.alignmentset.bam` files.


This unsupported method has the following steps:

## 0_splitRef.sh
*Split reference into two files*

Input: primary contigs and haplotigs from FALCON Unzip assembly.

Output: two contig sets of roughly equal size and length distribution.

### Usage: 
```bash
0_splitRef.sh myPrimaryContigs.fa myHaplotigs.fa
```

### Dependencies
[samtools](http://www.htslib.org/)

[falcon_tools](https://github.com/gconcepcion/falcon_tools) which is part of [FALCON](http://pb-falcon.readthedocs.io/en/latest/quick_start.html) installation

## 1_resequencing
User must perform two separate resquencing analyses using the contig sets generated in the previous step. Downstream steps depend on [SMRT Link](http://www.pacb.com/support/software-downloads/) job directories.

## 2_sharedReads.sh
*Generate list of reads mapped to both contig sets*

NOTE: This is the most time and compute intensive step. Samtools uses a lot of memory at this step so monitor your resource usage. Perhaps there is a better way to pull the read IDs?

Input: Two SMRT Link job paths, number of processors to use.

Output: Text file containing list of reads that are shared between the resulting BAM files of the resequencing jobs.

### Usage

```bash
2_sharedReads.sh SMRTLinkJobDir1 SMRTLinkJobDir2 12
```

### Dependencies

[samtools](http://www.htslib.org/)

[GNU parallel](https://www.gnu.org/software/parallel/)

## 3_merge_ReseqBAMs.sh
*Rewrite resulting BAM files from resequencing jobs to remove lower quality alignments by calling auxillary python script, "MergeBAMpair.py."*


NOTE: New BAM files are written to the local subdir directory, "newBAMs". 

### Dependencies

python 2.7

pysam, os, csv, argparse, collections, sys

[GNU parallel](https://www.gnu.org/software/parallel/)

