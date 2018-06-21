source /mnt/software/Modules/current/init/bash
module load parallel/20161222
#2016-12-31 version of gnu parallel
module load samtools

mkdir newBAMs

# remove worse alignment for shared reads
# rewrite BAM files
$nproc = $1
parallel --verbose --link -j $nproc MergeBAMpair.py {1} {2} :::: BAM1.fofn BAM2.fofn


### optional  next steps

# file of filename for new BAMs
#find ./newBAMs -name "*.bam" | sort > new_BAMs.fofn

# index new BAMs with samtools
#parallel --verbose -j $nproc samtools index {} :::: new_BAMs.fofn

# samtools merge
#samtools merge -f -b new_BAMs.fofn > merged.bam

# index new BAMs with pbindex
#parallel --verbose -j $nproc pbindex {} :::: new_BAMs.fofn

# new pb alignment dataset
#dataset create --type AlignmentSet --name merged_alignmentset merged_alignmentset.xml new_BAMs.fofn 

