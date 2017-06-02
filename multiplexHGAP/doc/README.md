# Multiplexed Analysis Explained
In general, running demultiplexed analyses in SMRT Analysis requires three basic steps: 

    1) Barcode
    2) Demultiplex
    3) Analysis

This document is intended to explain these three steps (**without** SMRT-Link) so that users can generate custom command-line workflows as necessary. These are essentially the steps that are used in the scripts found in this repository (via SMRT-Link), but here they are explained such that users can implement their own scripts if desired.

## Terminology
* **Barcoding**
    Label BAM rows in an UNALIGNED SubreadSet with tags 'bc' and 'bq'.
    * Input: subreads.bam, scraps.bam (or a subreadset.xml including both)
    * Output: subreads.bam, scraps.bam 
* **Demultiplex** Symbolic or physical separation of barcoded datasets by barcode values.
    The *dataset* abstraction allows symbolic separation using XML filter files.  Physical separation into separate BAM files per barcode is possible with a separate *consolidation* step to generate new BAMs per barcode.
* **Analysis**
    Any SMRT Analysis application can use the symbolic (XML) or physical (BAM plus XML) output of demultiplexing as in input for analysis. 

## Pbsmrtpipe templates
Users can easily see available applications using pbsmrtpipe with the following:

    pbsmrtpipe show-templates

To generate a preset template for a pipeline:

    pbsmrtpipe show-template-details -j hgap4_template.json pbsmrtpipe.pipelines.polished_falcon_fat

## Barcode
Barcoding in SMRT Analysis v3+ is accomplished with the tool `bam2bam`.  On the command line, this can be done with either `bam2bam` directly, or using `pbsmrtpipe` to take advantage of chunking and parallelization of the analysis.
    
Direct barcoding:

    bam2bam -j 8 -b 8 \
            --barcodes myBarcodes.barcodeset.xml \
            -o myBarcodedReads \
            --scoreMode symmetric \
            <movie>.subreadset.xml
    
    or (using bams):

    bam2bam -j 8 -b 8 \
            --barcodes myBarcodes.barcodeset.xml \
            -o myBarcodedReads \
            --scoreMode symmetric \
            <movie>.subreads.bam <movie>.scraps.bam

Using pbsmrtpipe:

    pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_barcode \
                           -o outputDirectory \
                           -e "eid_subread:<movie>.subreadset.xml" \
                           -e "eid_barcode:myBarcodes.barcodeset.xml"

To create a subreadset.xml from arbitrary bam(s) for use with PacBio applications:

    dataset create --type SubreadSet \
                   --name mySubreadsetName \
                   mySubreadSetName.subreadset.xml \
                   a.bam b.bam [...]

Note that for generating **new** barcode calls in the first step, **both** the subreads.bam and scraps.bam must be included in the subreadset.xml definition. 

PacBio bam files which are already barcoded (see bam header) do not need to include the scraps.bam file for downstream applications which use the barcode tags.  See [PacBioFileFormats](https://github.com/PacificBiosciences/PacBioFileFormats/blob/4.0/BAM.rst#bam-filename-conventions) for more information.

## Demultiplex
Generally, it is only necessary to generate a dataset xml file defining each subset corresponding to the barcodes in the original data, rather than physically separating the barcoded BAM files.  Each application will extract the requisite raw reads for analysis, given the dataset definition.

The `bam2bam` application in the previous step will generate a `subreadset.xml` file in the output location which can be passed to the demultiplexing step.  Alternately, `pbsmrtpipe` will generate a *gathered* `file.subreadset.xml` dataset in the directory `<outDir>/tasks/pbcoretools.tasks.gather_subreadset-1/` which can be passed to the next step. 

    dataset split --barcodes \
                  --outdir mySplitDatasets/ \
                  myBarcodedReads.subreadset.xml

The result of this step will be a set of `barcoded.chunk##.subreadset.xml` files in the output directory, each of which corresponds to a filter on the barcoded dataset for one barcode in the barcodeset that was used in the barcoding step above.  

Please note that the **chunk##** does **NOT** correspond to the barcode or barcode index.  The chunk numbers are named by the number of records/rows in each subset in descending order. The barcode *index pair* for each subreadset can be obtained by looking inside the XML file for the `bc=[#,#]` filter, or by calling the following command:

    dataset summarize barcoded.chunk#.subreadset.xml

### Filter
An additional filtering step is recommended for most applications to prevent inclusion of false-positive calls in analysis.  This step will filter barcode calls based on a minimum barcode score, here we recommend a minimum score of 30-40 for most applications.

    dataset filter barcoded.chunk#.subreadset.xml barcoded.chunk#.filtered.subreadset.xml "bq>=40"

## Analysis
The resultant (filtered) demultiplexed subreadsets can be used as entry subreadsets to any pbsmrtpipe application.  Generate and edit the preset json files according to steps above.

    pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.polished_falcon_fat \
                           -o myOutDir \
                           --preset-json myPresets.json \
                           -e "eid_subread:barcoded.chunk#.filtered.subreadset.xml"

    pbsmrtpipe pipeline-id pbsmrtpipe.pipelines.sa3_ds_resequencing_fat \
                           -o myOutDir \
                           --preset-json myPresets.json \
                           -e "eid_subread:barcoded.chunk#.filtered.subreadset.xml" \
                           -e "eid_ref_dataset:myReference.xml"
