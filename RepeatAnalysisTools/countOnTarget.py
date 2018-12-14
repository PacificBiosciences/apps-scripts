#! /home/UNIXHOME/jharting/anaconda2/bin/python

import pandas as pd
import sys,os,pysam
from resources.utils import readBED,getZmw

FLOATFORMAT='%.4f'

def main(parser):
    args = parser.parse_args()

    #get targets
    targets           = readBED(args.inBED,names=['ctg','start','end','name'])
    targets['length'] = targets.eval('end - start')
    targets.set_index('name',inplace=True)
    #open alignments
    bam      = pysam.AlignmentFile(args.inBAM,'rb')
    counter  = makeCounter(bam)
    #count alignments passing 'good' filter
    try:
        counts = pd.Series({ target : counter(row.ctg,row.start,row.end)
                            for target,row in targets.iterrows()})\
                   .rename('onTargetZMWs')
    except ValueError,e:
        #try to catch input bams without .bai
        raise countOnTarget_Exception(e)
    #calculate enrichment and build result table
    nZMWs,avgReadLength   = getReadStats(bam)
    genomeSize            = sum(bam.lengths)
    targets['expected']   = (nZMWs * (avgReadLength - targets.length) / genomeSize).clip_lower(0)
    results               = targets.join(counts)
    results['enrichment'] = results.eval('onTargetZMWs / expected')
    #write output
    results.to_csv('{}/onTargetCounts.tsv'.format(args.outdir),
                   float_format=FLOATFORMAT,header=True)

    return results

def isGoodAlignment(rec,start,stop):
    #spanning reads.  
    #Add other filters as needed for map quality, etc here
    return rec.reference_start < start and rec.reference_end > stop

def makeCounter(bam,**kwargs):
    def counter(ctg,start,stop):
        return len({getZmw(rec.query_name) 
                    for rec in bam.fetch(ctg,start,stop) 
                    if isGoodAlignment(rec,start,stop,**kwargs)})
    return counter

def getReadStats(bam):
    n = 0
    mean = 0.0
    bam.reset()
    for rec in bam:
        n += 1
        mean += (rec.query_length - mean) / n
    return n,mean

class countOnTarget_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='countOnTarget.py', description='Generate table of ZMW counts per target')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='BAM file of aligned reads.  Must have .bai index')
    parser.add_argument('inBED', metavar='inBED', type=str,
                    help='BED file with targets')
    parser.add_argument('-o,--outdir', dest='outdir', type=str, default=os.getcwd(),
                    help='directory to save output file.  default cwd.')
    
    try:
        result = main(parser)
    except countOnTarget_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
