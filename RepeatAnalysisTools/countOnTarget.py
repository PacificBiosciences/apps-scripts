#! /home/UNIXHOME/jharting/anaconda2/bin/python

import pandas as pd
import sys,os,pysam

def main(parser):
    args = parser.parse_args()

    #get targets
    targets  = readBED(args.roiBED)
    #open alignments
    bam      = pysam.AlignmentFile(args.inBam,'rb')
    counter  = makeCounter(bam)
    #count alignments passing 'good' filter
    try:
        counts = pd.Series({ row["name"] : counter(row.ctg,row.start,row.end)
                            for i,row in targets.iterrows()})\
                   .rename('onTargetZMWs')
    except ValueError,e:
        #try to catch input bams without .bai
        raise countOnTarget_Exception(e)

    counts.to_csv('{}/onTargetCounts.tsv'.format(args.outdir),
                  sep='\t',header=True)

    return counts


def readBED(bedfile):
    return pd.read_table(bedfile,sep='\s+',names=['ctg','start','end','name'])

def isGoodAlignment(rec,start,stop):
    #spanning reads.  
    #Add other filters as needed for map quality, etc here
    return rec.reference_start < start and rec.reference_end > stop

def makeCounter(bam):
    def counter(ctg,start,stop):
        return len({getZmw(rec.query_name) for rec in bam.fetch(ctg,start,stop) if isGoodAlignment(rec,start,stop)})
    return counter

def getZmw(readname):
    return '/'.join(readname.split('/')[:2])

class countOnTarget_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='countOnTarget.py', description='Generate table of ZMW counts per target')
    parser.add_argument('inBam', metavar='inBam', type=str,
                    help='BAM file of aligned reads.  Must have .bai index')
    parser.add_argument('roiBED', metavar='roiBED', type=str,
                    help='BED file with targets')
    parser.add_argument('-o,--outdir', dest='outdir', type=str, default=os.getcwd(),
                    help='directory to save output file.  default cwd.')
    
    try:
        result = main(parser)
    except countOnTarget_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
