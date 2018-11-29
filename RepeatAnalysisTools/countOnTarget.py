#! /home/UNIXHOME/jharting/anaconda2/bin/python

import pandas as pd
import numpy as np
import os
from pbcore.io import IndexedBamReader


def main(parser):
    args = parser.parse_args()

    #get targets
    targets  = readBED(args.roiBED)
    #open alignments
    bam      = IndexedBamReader(args.inBam)
    #read alignment pos from index
    columns  = ['tId','tStart','tEnd']
    alnTable = pd.DataFrame(np.array([getattr(bam.index,col) 
                                      for col in columns]).T,
                            columns=columns)
    #asign names to index ids
    refMap          = dict(zip(bam.referenceInfoTable.ID,bam.referenceInfoTable.Name))
    alnTable['ctg'] = alnTable.tId.map(refMap)

    #generate results
    qry = 'ctg==@row.ctg & tStart<@row.start & tEnd>@row.end'
    counts = pd.Series({row["name"]:len(alnTable.query(qry))
                        for i,row in targets.iterrows()})\
               .rename('onTargetZMWs')

    counts.to_csv('{}/onTargetCounts.tsv'.format(args.outdir),
                  sep='\t',header=True)

    return counts


def readBED(bedfile):
    return pd.read_table(bedfile,sep='\s+',names=['ctg','start','end','name'])

class countOnTarget_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='countOnTarget.py', description='Generate table of ZMW counts per target')
    parser.add_argument('inBam', metavar='inBam', type=str,
                    help='BAM file of aligned reads.  Must have .pbi index')
    parser.add_argument('roiBED', metavar='roiBED', type=str,
                    help='BED file with targets')
    parser.add_argument('-o,--outdir', dest='outdir', type=str, default=os.getcwd(),
                    help='directory to save output file.  default cwd.')
    
    try:
        result = main(parser)
    except countOnTarget_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
