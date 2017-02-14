#! /bin/bash

NAME=$1                               #Base name for jobs/subreadsets
BARCODESUBREADSET=$(readlink -f $2)   #Barcoded subreadset xml
BARCODEFASTA=$(readlink -f $3)        #Barcode fasta used for barcoding subreadset.  Used for naming.
GENOMESIZE=$4                         #Estimated genome size in bp
OUTDIR=$(readlink -f $5)              #Directoy to place outputs
HOST=$6                               #SMRT Link host name
PORT=$7                               #SMRT Link service port 

MINBCSCORE=45    #Minimum barcode SW alignment score [0-100].  Recommended 45.
MINSUBREADS=5000 #Minimum number of records per barcode (else do not submit job)
SLEEP=5          #minutes to sleep between submitting jobs

PYTHON=/path/to/smrtlink/installdir/private/otherbins/all/bin/python
SPLITPROG=/path/to/splitBarcodeUpload.py
JOBPROG=/path/to/multiplexSubmit.py

$PYTHON $SPLITPROG $BARCODESUBREADSET -n $NAME         \
                                      -b $BARCODEFASTA \
                                      -o $OUTDIR       \
                                      -f $MINBCSCORE   \
                                      -r $MINSUBREADS  \
                                      --host $HOST     \
                                      --port $PORT

$PYTHON $JOBPROG -o $OUTDIR                              \
                 -S $OUTDIR/uploaded_subreadsets.csv     \
                 --HGAP_GenomeLength_str $GENOMESIZE     \
                 -w $SLEEP                               \
                 --host $HOST                            \
                 --port $PORT 
