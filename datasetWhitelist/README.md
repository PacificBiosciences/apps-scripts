# Whitelisting Subread Datasets for SMRT Analysis v4.0 through v5.0.1 
## *DEPRECATED for SMRT analysis v5.1* 

Starting with SMRT Analysis v5.1, whitelisting can be done natively in pbcore with the `dataset` tool by passing a file of read names or ZMWs as the filter argument.  Set intersetions are handled by pbcore 'under the hood'

For subreads:

    dataset filter in.dataset.xml filtered.dataset.xml 'qname=/path/to/textFile1.txt'

Where `textFile1.txt` is a simple text file of PacBio subread names `[movie]/[holenumber]/[qStart]_[qEnd]`, one per line

For ZMWs:

    dataset filter in.dataset.xml filtered.dataset.xml 'zm=/path/to/textFile2.txt'

Where `textFile2.txt` is a simple text file with sequencing hole numbers, one per line (Note that this assumes a single smrtcell/movie).
This method will create a new dataset (view) of the filtered data without copying raw sequence data into a new BAM(s).

To blacklist reads/holes, use the `!=` operator in the above filter commands.


Alternatively, the tool `bamsieve` in the SMRT Analysis v5.1 installation can be used to create new filtered BAM files from white-/blacklists:

    bamsieve --whitelist textFile3.txt in.bam out.filtered.bam

See `bamsieve --help` for more information.

Filtered subreadsets can be uploaded to the SMRT-Link server using `pbservice import-dataset` 

## SMRT Analysis v4.0-5.0
Users wishing to specify exactly which ZMWs or subreads are used in SMRT Analysis applications can apply "whitelist" filtering to datasets using the script here.  Whitelisted dataset xml files can be used as input to any SMRT Analysis application which accepts a dataset as input, including SMRT-Link protocols or command-line tools.
## Dependencies
[pbcore](https://github.com/PacificBiosciences/pbcore)

## Whitelist ZMWs

    SSET=/path/to/mydata.subreadset.xml

Include all subreads from any ZMW identified by reads in `myReads.fasta`

    python datasetWhitelist.py $SSET myReads.fasta              \
                               -o myFilteredData.subreadset.xml \
                               -n myNewFilteredSubreads
    
Include all subreads from any ZMW identified by BLASR mapping

    blasr myReads.fasta myRef.fasta | python datasetWhitelist.py -l $SSET -o myFilteredData.subreadset.xml 

## Whitelist Subreads
Include only subreads explicitly listed in `myReads.fasta`

    python datasetWhitelist.py -s $SSET myReads.fasta           \
                               -o myFilteredData.subreadset.xml \
                               -n myNewFilteredSubreads

Include only subreads explicitly identified by BLASR mapping

    blasr myReads.fasta myRef.fasta | python datasetWhitelist.py -s -l $SSET -o myFilteredData.subreadset.xml

## Blacklist
The inverse of the above opperations can be done by adding the -i option to 'blacklist' the input read names. 

    python datasetWhitelist.py -s -i $SSET myReads.fasta        \
                               -o myFilteredData.subreadset.xml \
                               -n myNewFilteredSubreads

## Check the --help

    python datasetWhitelist.py -h

    usage: datasetWhitelist.py [-h] -o,--outXml OUTXML [-n,--name NAME]
                               [-s,--subreads] [-l,--list] [--noUuid]
                               [-i,--inverted]
                               inXml [inFile]
    
    Generate a whitelisted dataset xml from an input fasta of reads (subreads or
    ccs) or file of readnames
    
    positional arguments:
      inXml               input pacbio subread dataset.
      inFile              file with read names (e.g. fasta,blasr output,text file
                          of names). default stdin
    
    optional arguments:
      -h, --help          show this help message and exit
      -o,--outXml OUTXML  output xml.
      -n,--name NAME      name for filtered dataset. default keep original name
      -s,--subreads       whitelist subread names instead of zmws (will only work
                          with subread names, not ccs names). default false.
      -l,--list           input names as text list. default false.
      --noUuid            do not generate a new uuid for the dataset. default true
                          (create new uuid).
      -i,--inverted       invert list. i.e. 'Blacklist' the input values from the
                          dataset. default false.
