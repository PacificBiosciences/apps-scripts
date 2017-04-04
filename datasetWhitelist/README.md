# Whitelisting Subread Datasets for SMRT Analysis v4.0+
Users wishing to specify exactly which ZMWs or subreads are used in SMRT Analysis applications can apply "whitelist" filtering to datasets using the script here.  Whitelisted dataset xml files can be used as input to any SMRT Analysis application which accepts a dataset as input, including SMRT-Link protocols or command-line tools.

Filtered subreadsets can be uploaded to the SMRT-Link server using `pbservice import-dataset` 

## Dependencies
[pbcore](https://github.com/PacificBiosciences/pbcore)

## Whitelist ZMWs

    SSET=/path/to/mydata.subreadset.xml

Include all subreads from any ZMW identified by reads in `myReads.fasta`

    python datasetWhitelist.py $SSET myReads.fasta -o myFilteredData.subreadset.xml
    
Include all subreads from any ZMW identified by BLASR mapping

    blasr myReads.fasta myRef.fasta | python datasetWhitelist.py --list $SSET -o myFilteredData.subreadset.xml 

## Whitelist Subreads
Include only subreads explicitly listed in `myReads.fasta`

    python datasetWhitelist.py --subreads $SSET myReads.fasta -o myFilteredData.subreadset.xml

Include only subreads explicitly identified by BLASR mapping

    blasr myReads.fasta myRef.fasta | python datasetWhitelist.py --subreads --list $SSET -o myFilteredData.subreadset.xml

## Check the --help

    python datasetWhitelist.py -h

    usage: datasetWhitelist.py [-h] -o,--outXml OUTXML [-s,--subreads] [-l,--list]
                               inXml [inFile]
    
    Generate a whitelisted dataset xml from an input fasta of reads (subreads or
    ccs) or file of readnames (e.g. blasr)
    
    positional arguments:
      inXml               input pacbio subread dataset.
      inFile              file with read names (e.g. fasta,blasr output,text file
                          of names). default stdin
    
    optional arguments:
      -h, --help          show this help message and exit
      -o,--outXml OUTXML  output xml.
      -s,--subreads       whitelist subread names instead of zmws (will only work
                          with subread names, not ccs names). default false.
      -l,--list           input names as text list. default false.
