# Reseq_BigGenome
## Resequencing (polishing) hack for draft genomes larger than 4.29 Gb

The resequencing pipeline is recommended for polishing genomes assembled with PacBio data to increase base accuracy. Resequencing uses [BLASR](https://github.com/PacificBiosciences/blasr) to map raw reads to the draft reference and arrow for [genomic consensus](https://github.com/PacificBiosciences/GenomicConsensus) and can be accessed with the [SMRT Link GUI](http://www.pacb.com/wp-content/uploads/SMRT_Link_User_Guide.pdf) or [pbsmrtpipe](http://pbsmrtpipe.readthedocs.io/en/master/getting_started.html#basic-resequencing) at the command line. BLASR cannot index reference sequences longer than 2^32 (4.29 Gb).

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

This is the most time and compute intensive step. Samtools uses a lot of memory at this step so monitor your resource usage. Perhaps there is a better way to pull the read IDs?

Input: Two SMRT Link job paths, number of processors to use.

Output: Text file containing list of reads that are shared between the resulting BAM files of the resequencing jobs.

### Usage

```bash
2_sharedReads.sh SMRTLinkJobDir1 SMRTLinkJobDir2 12
```

### Dependencies

[samtools](http://www.htslib.org/)

[GNU parallel](https://www.gnu.org/software/parallel/)

## 3_mergeReseqJobs.sh
*Modify resulting BAM files from Resequencing Jobs to remove lower quality alignments*

For each of these reads, the BAM files are queried and the alignment with the better (more negative) BLASR score is retained while the other is omitted. NOTE: New BAM files are written to the local directory.


