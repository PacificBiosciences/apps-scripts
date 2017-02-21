# Multiplexed Microbial Assembly with SMRT-Link
Running multiplexed samples for assembly directly in the SMRT-Link GUI is not currently supported as an automated process.  The tools to do so are readily available both in the GUI and on the command line.  These instructions will work with **SMRT Analysis v4.0** and later versions.

## Target Audience
This wiki is intended for users who have multiplexed datasets to be passed through SMRT-Link applications on a per-barcode basis.  The most common example is multiplexed assembly of small (microbial) genomes using the HGAP 4 application.  Other downstream applications, e.g. resequencing per barcode, can also be run using the instructions here.  

Note that LAA2 and CCS2 applications handle barcoding within the workflow of a single job and users **do not** need this method for those analyses.

## Overview
This tutorial describes two methods for running SMRT-Link applications on a per-barcode basis.  
* Manual Setup (GUI)
* Automated Submission (Command Line)

### Barcoded SubreadSet
Both methods require running the standard Barcoding application prior to the steps listed here.  The output barcoded subreadset from the Barcoding application is the input for this tutorial.

### Collating Results
A final option to collate reports from multiple SMRT-Link jobs (using the same application) is provided at the end for command-line users.

## Manual HGAP per Barcode Using GUI
Step 1: Select barcoded data. From the SMRT-Link home page, 
    
    SMRT Analysis > Create New Analysis > Select your barcoded dataset > Next

Note that subreadsets with barcodes have the label "(barcoded)" added to the "Name" field.

Step 2: Select ``Assembly (HGAP4)`` from the ``Analysis Application`` drop-down menu and set the approximate ``Genome Length``.

Step 3: In ``Advanced Analysis Parameters > Filters to add to the DataSet``, add comma-separated filters to the input dataset

    Example A (symmetric):
    rq >= 0.7, bcf = 0, bcr = 0, bq >= 45

    Example B (asymmetric,universal):
    rq >= 0.7, bcf = 0, bcr = 1, bq >= 45

* **rq** : Read Quality filter.  Default setting is usually sufficient.
* **bcf**: Barcode Forward.  The 0-indexed position of the forward barcode in the BarcodeSet used for barcoding*
* **bcr**: Barcode Reverse.  The 0-indexed position of the reverse barcode in the BarcodeSet used for barcoding
* **bq** : Barcode Quality Score.  SW alignment score of the barcode on integer scale of 0-100.  Recommended minimum 45.

![Advanced Analysis Parameters Screenshot](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/images/manual_filtering_screenshot.png)

Step 4: Change any other parameters as desired.  Click **Start**.

Step 5: Manually Repeat steps 1-4 for each barcode in the dataset.

### Indexing Barcodes
You will need to take a look at the fasta file underlying the BarcodeSet that you used for barcoding in order to identify indices with the barcode names.  For the two example filter parameters above and the following snippet of the first four records in the barcode fasta,

    >lbc1
    TCAGACGATGCGTCAT
    >lbc2
    CTATACATGACTCTGC
    >lbc3
    TACTAGAGTAGCACTC
    >lbc4
    TGTGTATCAGTACATG

the examples correspond to the following barcodes:

``Example A (symmetric):``**``lbc1``**``(index 0 on both sides of the insert)``

``Example B (asymmetric,universal):``**``lbc1--lbc2``**``(index 0 with the forward primer, index 1 with the reverse primer)``


## Auotomated HGAP Job Submission per Barcode using Python

Automating barcoded job submission can be done using the tools and Python API in the SMRT-Link installation.  The scripts needed for this tutorial can be found below.  The scripts use the [pbcommand](https://github.com/PacificBiosciences/pbcommand) library and the [pbservice](https://github.com/PacificBiosciences/SMRT-Link/wiki/pbservice) tool to interact with SMRT-Link services.  You will need to know the location of your SMRT-Link install to modify the scripts for your system. 

* [Bash quickstart (multiplexHGAP4.sh)](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/scripts/multiplexHGAP4.sh)
* [Splitting and Importing (splitBarcodeUpload.py)](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/scripts/splitBarcodeUpload.py)
* [Submitting Jobs (multiplexSubmit.py)](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/scripts/multiplexSubmit.py)
* [HGAP presets json (presets_template.json)](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/templates/presets_template.json)

### Quickstart to Automation

Step 1:  Modify ``multiplexSubmit.py`` to correspond to where you saved the ``presets_template.json``, your SMRT Analysis install directory, and the ``host`` and ``port`` (See 'About' link in the SMRT-Link browser Menu).

    #Set the following values at the top of the script
    PRESETS_TEMPLATE='/path/to/presets_template.json'
    PBSERVICE       ='/path/to/smrtlink/installdir/smrtcmds/bin/pbservice'
    DEFAULTHOST     ='smrtlink'
    DEFAULTPORT     =8081

Step 2: Modify ``multiplexHGAP4.sh`` to point to the SMRT-Link python interpreter, and the two scripts you downloaded above.

    #Set the following values in the script
    PYTHON=/path/to/smrtlink/installdir/private/otherbins/all/bin/python
    SPLITPROG=/path/to/splitBarcodeUpload.py
    JOBPROG=/path/to/multiplexSubmit.py

You might also consider setting the ``HOST`` and ``PORT`` values rather than passing them in as arguments to the script.

Step 3:  Run the script using a barcoded subreadset xml as input

    NAME='My Multiplexed Jobs'                         
    BARCODESUBREADSET=/path/to/barcoded.subreadset.xml  #Can use xml downloaded from SMRT-Link 
    BARCODEFASTA=/path/to/barcode.fasta                 #Barcode fasta used to score subreadset. Used for naming only.
    GENOMESIZE=4600000                                  #Estimated
    OUTDIR=/my/output/directory
    HOST=smrtlink                                       #See SMRT-Link 'About' page in browser
    PORT=8081

    bash multiplexHGAP.sh $NAME $BARCODESUBREADSET $BARCODEFASTA $GENOMESIZE $OUTDIR $HOST $PORT

The bash script runs the two python scripts. First run ``splitBarcodeUpload`` to generate subreadset.xml files, one for each barcode and import them to the SMRT-Link server.  Second run ``multiplexSubmit.py`` to start HGAP4 jobs, one for each imported barcode subreadset using default parameters and specified estimated genomesize.


After the script finishes, the output directory contains the following:

    ls /my/output/directory

    [NAME]_lbc*--lbc*.subreadset.xml            #Filtering xmls, one for each barcode
    [NAME]_barcodelbc*--lbc*_HGAP_presets.json  #HGAP presets,  one for each barcode
    uploaded_subreadsets.csv                    #Table of imported datasets.  (subsetName,subsetId,jobState)
    started_jobs.csv                            #Table of started jobs. (host,jobId,jobName,jobPath)

Jobs will usually still be running after the submit script finishes, and they will all run concurrently by default.  If you are concerned for overloading your cluster environment with too many jobs, you can set the ``SLEEP`` parameter in ``multiplexHGAP.sh`` to wait N minutes between submitting jobs.

### Splitting and Importing each barcode
The script ``splitBarcodeUpload.py`` can be used independently to upload barcoded subsets to you SMRT-Link server for further analysis in any downstream SMRT-Link application.  See ``python splitBarcodeUpload.py --help`` for more information.

Note that this script does *not* physically demultiplex the original subreads.bam.  Each imported subset is an xml file pointing to the original barcoded subreadset with an associated set of filters for barcode value and quality.  

See [documentation](http://pacbiofileformats.readthedocs.io/en/3.0/DataSet.html) for more information on the dataset model. 

### Multiple Job Sumbission
The script ``multiplexSubmit.py`` can also be used independently to submit multiple SMRT-Link jobs given a presets.json template and a csv file of input subreadsets IDs from your server.  Just change the template in the script to use this with other SMRT-Link applications.  The script will pick up the available options from the json file.

    PRESETS_TEMPLATE='/path/to/presets_template.json'

See ``multiplexSubmit.py --help`` after setting the preset template.

The json file linked to on this wiki is for HGAP4 and can be generated using ``pbsmrtpipe`` from the SMRT Analysis environment.  Generating a json file for resequencing can be done as follows:

    pbsmrtpipe show-template-details pbsmrtpipe.pipelines.sa3_ds_resequencing -j presets.json


## Comparing Results
The following script is provided as a way to collate results reports for multiple SMRT-Link jobs.  Note that this tool depends on the python module [pandas](http://pandas.pydata.org/) (not incuded in the SMRT Analysis python) as well as [pbcommand](https://github.com/PacificBiosciences/pbcommand).  

* [collateReports.py](https://github.com/PacificBiosciences/apps-scripts/blob/master/multiplexHGAP/scripts/collateReports.py)

This takes a csv file with the minimum columns ``host``, ``jobId``, and ``jobName`` and attempts to collate reports from the listed jobs for export to csv and excel formats.

    head started_jobs.csv

    host,jobId,jobName,jobPath
    smrtlink,1934,12plex_bsub (barcode=lbc34--lbc34)_HGAP,/path/to/jobs/001934
    smrtlink,1935,12plex_bsub (barcode=lbc40--lbc40)_HGAP,/path/to/jobs/001935
    smrtlink,1936,12plex_bsub (barcode=lbc29--lbc29)_HGAP,/path/to/jobs/001936
    smrtlink,1937,12plex_bsub (barcode=lbc38--lbc38)_HGAP,/path/to/jobs/001937
    smrtlink,1938,12plex_bsub (barcode=lbc70--lbc70)_HGAP,/path/to/jobs/001938
    
See ``collateReports.py --help`` for more information.

