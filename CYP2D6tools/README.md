# Tools for star-typing analysis of CYP2D6

Under Construction!

### Dependencies
 - [pysam](https://github.com/pysam-developers/pysam)
 - [mappy](https://pypi.org/project/mappy/)
 - [pandas](https://pandas.pydata.org/)
 - [numpy](https://numpy.org/)
 - [sqlalchemy](https://www.sqlalchemy.org/)


    usage: pbCYP2D6typer.py [-h] -r RUNNAME [-s SAMPLEMAP] [-i] [-f MINFRAC]
                            [-P {splice,map-pb,gaplenient,gapstrict}]
                            [--hifiSupport HIFISUPPORT] [--read_info READ_INFO]
                            [-p PREFIX] [--writeCSV] [--vcf] [--dbStore]
                            [consensusFasta [consensusFasta ...]]
    
    Call CYP2D6 star alleles from pbaa consensus
    
    positional arguments:
      consensusFasta        Fasta file(s) of pbAA consensus outputs
    
    optional arguments:
      -h, --help            show this help message and exit
    
    Input Options:
      -r RUNNAME, --runName RUNNAME
                            Sequencing run ID. Required
      -s SAMPLEMAP, --sampleMap SAMPLEMAP
                            CSV file mapping biosample name to barcode. Must have
                            fields ['Barcode', 'Bio Sample Name']. Default None
      -i, --ignoreMissing   Do not report empty rows for missing samples from
                            sample map. Default False
      -f MINFRAC, --minFrac MINFRAC
                            Ignore failed clusters below minFrac. Default 0.01
      -P {splice,map-pb,gaplenient,gapstrict}, --preset {splice,map-pb,gaplenient,gapstrict}
                            Alignment preset for mappy aligner. Choose "splice"
                            for expected large deletions. Default gapstrict
      --hifiSupport HIFISUPPORT
                            Add per-variant read support depth. Path to hifi fastq
                            used for clustering. Requires read_info option to be
                            set. Default None
      --read_info READ_INFO
                            Path to pbAA read_info table for the run. Required for
                            hifi-support option. Default None
    
    Output Options:
      -p PREFIX, --prefix PREFIX
                            Output prefix. Default ./cyp2d6typer
      --writeCSV            Write intermediate CSV files. Default False
      --vcf                 Output VCF file. HiFi support and read_info required.
                            Default False
      --dbStore             Maintain sample records in local db. Default False

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.

FOR RESEARCH USE ONLY. NOT FOR USE IN DIAGNOSTICS PROCEDURES.
