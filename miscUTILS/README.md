# Description
This repo contains miscellaneous stand-alone scripts for working with PacBio Data

## getJobData.py
Access downloadable data from SMRT Link jobs via the command line.  Can be used to print SMRT Link job output data names and locations, or to copy the outputs to a new location.
### Dependencies
 - pbcommand
### Usage
    $ python getJobData.py -h
    usage: getJobData.py [-h] [-o,--outDir OUTDIR] [-n,--names] [--host HOST]
                         [--port PORT]
                         jobNumber [dataName [dataName ...]]
    
    Get job outputs from smrtlink datastore
    
    positional arguments:
      jobNumber           SMRT Link job number
      dataName            name of datasetset to get. default None. (use -n to see
                          available names)
    
    optional arguments:
      -h, --help          show this help message and exit
      -o,--outDir OUTDIR  output destination directory. default cwd
      -n,--names          print names of available items in datastore. default
                          False
      --host HOST         smrtlink host name. default smrtlink
      --port PORT         smrtlink services port. default 8081
### Examples
    Print description and paths to output files

    $ python getJobData.py --host smrtlink-release \
                           --port 9091 \
                           -n \
                          12258

    Amplicon Consensus Report               /path/to/012258/tasks/pbreports.tasks.amplicon_analysis_consensus-0/consensus_report.json
    Consensus Sequence Statistics           /path/to/012258/tasks/pbcoretools.tasks.gather_csv-1/file.csv
    Chimeric/Noise Sequences by barcode     /path/to/012258/tasks/pbcoretools.tasks.split_laa_fastq-0/chimera_fastq.gz
    Consensus Amplicons                     /path/to/012258/tasks/pbcoretools.tasks.split_laa_fastq-0/consensus_fastq.gz
    Input Molecule Report CSV               /path/to/012258/tasks/pbcoretools.tasks.gather_csv-2/file.csv
    Consensus Amplicons                     /path/to/012258/tasks/pbcoretools.tasks.gather_fastq-1/file.fastq
    LAA Input Report                        /path/to/012258/tasks/pbreports.tasks.amplicon_analysis_input-0/amplicon_input_report.json
    Chimeric/Noise Sequences                /path/to/012258/tasks/pbcoretools.tasks.gather_fastq-2/file.fastq
    Master Log                              /path/to/012258/logs/master.log
    Analysis Log                            /path/to/012258/logs/pbsmrtpipe.log


    Copy outputs to local directory

    $ python getJobData.py --host smrtlink-release \
                           --port 9091 \
                           12258 \
                           'Amplicon Consensus Report' 'Input Molecule Report CSV' \
                           -o /my/output/directory

    'Amplicon Consensus Report'             =>      /my/output/directory/consensus_report.json
    'Input Molecule Report CSV'             =>      /my/output/directory/file.csv

## pbiPrintColumn.py
Print out values from .pbi index to stdout.  Useful for generating quick statistics for PacBio BAM data.
### Dependencies
 - pbcore
### Example
    $ python pbiPrintColumn.py demux.bam.pbi holeNumber bcQual qStart qEnd | head -10
    108989,46,13993,14536
    108989,46,14618,18457
    108989,46,18629,21091
    108989,46,21168,23475
    108989,46,23561,24022
    109019,70,1210,2891
    109019,70,2971,5223
    109019,70,5298,7568
    109019,70,7644,9846
    109019,70,9923,11915

Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
