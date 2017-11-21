Reseq_BigGenome
Resequencing (Polishing) Hack for Genomes larger than 4.29 Gb

The resequencing pipeline is recommended for polishing genomes assembled with PacBio data to increase base accuracy. Resequencing uses [BLASR](https://github.com/PacificBiosciences/blasr) to map raw reads to the draft reference and arrow for [genomic consensus](https://github.com/PacificBiosciences/GenomicConsensus) and can be accessed with the [SMRT Link GUI](http://www.pacb.com/wp-content/uploads/SMRT_Link_User_Guide.pdf) or [pbsmrtpipe](http://pbsmrtpipe.readthedocs.io/en/master/getting_started.html#basic-resequencing) at the command line. BLASR cannot index reference sequences longer than 2^32 (4.29 Gb).

This unsupported method has the following steps:

## 0_splitRef.sh
*Split reference into two files*

Input: primary contigs and haplotigs from FALCON Unzip assembly.
Output: two contig sets of roughly equal size and length distribution.

### usage: 
```bash
0_splitRef.sh myPrimaryContigs.fa myHaplotigs.fa
```

## 1_resequencing
User must perform two separate resquencing analyses using the contig sets generated in the previous step.

## 2_sharedReads.sh
*Generate list of reads mapped to both contig sets*

Input: Two SMRT Link job paths and number of processors to use
Output: Text file containing list of reads that are shared between the resulting BAM files of the resequencing jobs.

```bash
2_sharedReads.sh /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/userdata/jobs_root/076/076209 \
                  /pbi/dept/secondary/siv/smrtlink/smrtlink-beta/smrtsuite_166987/userdata/jobs_root/076/076210 \
                  24
```

## 3_mergeReseqJobs.sh
*Modify resulting BAM files from Resequencing Jobs to remove lower quality alignments*

For each of these reads, the BAM files are queried and the alignment with the better (more negative) BLASR score is retained while the other is omitted. NOTE: New BAM files are written to the local directory.


