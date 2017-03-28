# Whitelisting Subread Datasets for SMRT Analysis v4.0+
Users wishing to specify exactly which ZMWs or subreads are used in SMRT Analysis applications can apply "whitelist" filtering to datasets using the script here.  Whitelisted dataset xml files can be used as input to any SMRT Analysis application which accepts a dataset as input, including SMRT-Link protocols or command-line tools.

## Dependencies
[pbcore](https://github.com/PacificBiosciences/pbcore)

## Whitelist ZMWs
Include all subreads from any ZMW identified by reads in `myReads.fasta`
    python datasetWhitelist.py myData.subreadset.xml myReads.fasta myFilteredData.subreadset.xml 

## Whitelist Subreads
Include only subreads explicitly listed in `myReads.fasta`
    python datasetWhitelist.py --subreads myData.subreadset.xml myReads.fasta myFilteredData.subreadset.xml
