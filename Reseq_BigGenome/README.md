# PacBio_Reseq_BigGenome
Resequencing (Polishing) Hack for Genomes larger than 4.29 Gb

The resequencing pipeline is recommended for polishing genomes assembled with PacBio data. Resequencing uses BLASR to map raw reads to the draft reference and arrow for genomic consesnsus and can be accessed with the SMRT Link GUI or pbsmrtpipe at the command line. BLASR can not index reference sequences longer than 2^32 (4.29 Gb).

This unsupported method has three components:


