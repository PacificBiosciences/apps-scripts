#! /usr/bin/python

import numpy as np
import pandas as pd
import re,os,sys
import mappy as mp
import pysam
from tempfile import NamedTemporaryFile
from resources.utils import extractRepeat,\
                            getPositions,\
                            countAlignments,\
                            readBED,\
                            writeFasta
from resources.plotting import waterfallPlot,\
                               countPlot2

ALIGNFILTER=0x900
#Plotting resolution
DPI=400
#column name ofr target insert
SIZECOL='Insert Size (bp)'

def main(parser):
    args = parser.parse_args()
    
    #Read inputs
    bam  = pysam.AlignmentFile(args.inBAM)
    try:
        bed  = readBED(args.inBED,names=['ctg','start','end','name','motifs'])
    except Exception as e:
        #catch missing columns in bed
        raise RepeatAnalysisReport_Exception(e.message)
    ref  = pysam.FastaFile(args.reference)

    #function to make output file names
    s           = args.sample + '.' if args.sample else ''
    outfileName = lambda name,ext: '{d}/{s}{n}.{e}'.format(d=args.outDir,s=s,n=name,e=ext)
    
    targets,insertLength,summaries = {},{},{}
    for i,row in bed.iterrows():
        print 'Processing repeat regions: {n} -- {m}'.format(n=row['name'],m=row.motifs) 
        #cut out flanks from reference and use to identify repeat region
        #keep track of tmp reference file to remove it later
        #adjust start/stop for 0-based pysam coords
        aligner,tmp   = getFlankAligner(ref,row.ctg,row.start-1,row.end,
                                        flanksize=args.flanksize)
        #extract sequence between flanks
        repeatRegions = pd.DataFrame([(rec.query_name,) + extractRepeat(rec.query_sequence,aligner)
                                      for rec in bam.fetch(row.ctg,row.start,row.end)
                                      if not (rec.flag & ALIGNFILTER)],
                                      columns=['read','start','stop','subsequence'])
        repeatRegions.read = repeatRegions.read + '/' \
                                                + repeatRegions.start.astype(str) \
                                                + '_' \
                                                + repeatRegions.stop.astype(str) 
        repeatRegions[SIZECOL] = repeatRegions.subsequence.str.len()
        #make empty df for targets with no hits
        if not len(repeatRegions):
            repeatRegions    = pd.DataFrame(columns=['read','subsequence',SIZECOL])
        #sort records by length of repeat sequence
        repeatRegions.sort_values(SIZECOL,ascending=False,inplace=True)
        repeatRegions.reset_index(drop=True,inplace=True)
        #count alignments and filter for spanning reads
        summary,filtered = countAlignments(repeatRegions)
        #count repeat motifs
        motifs = row.motifs.split(',')
        try:
            #get repeat positions in sequence for each motif
            motifDfs = [pd.concat(filtered.set_index('read',append=True)\
                                          .subsequence.map(getPositions(motif))\
                                          .to_dict())\
                          .reset_index(level=2,drop=True)\
                          .reset_index()\
                          .rename(columns={'level_0':'idx','level_1':'readName'})
                        for motif in motifs]
            df = pd.concat(motifDfs,ignore_index=True)
        except ValueError:
            #catch targets with no motifs found, set df to empty w/columns
            df = pd.DataFrame(columns=['idx','readName','position','motif'])
        #add results to containers for concating later
        targets[row['name']]   = df
        summaries[row['name']] = summary
        insertLength.update(dict(filtered[['read',SIZECOL]].values))
        #write extracted regions
        if len(filtered):
            writeFasta(outfileName('exractedSequence_%s'%row['name'],'fasta'),
                       filtered[['read','subsequence']].values)

    #remove temp flank ref
        os.remove(tmp.name)
    #concat results
    repeatDf  = pd.concat(targets,names=['target','foo'])\
                  .reset_index('foo',drop=True)
    summaryDf = pd.concat(summaries).unstack()
    summaryDf.index.name = 'target'
    
    print 'Plotting Figures'
    #first motif for each tareet is primary
    primaryMotifs = zip(bed['name'],bed.motifs.str.split(',',expand=True)[0].values)
    for target,motif in primaryMotifs:
        try:
            data   = repeatDf.loc[target]
            #waterfall
            g = waterfallPlot(data) #target=target  (for label)
            g.savefig(outfileName('.'.join([target,'waterfall']),'png'),dpi=DPI)
            #histograms
            #plot main repeat motif
            counts = data.query('motif==@motif')\
                         .groupby('idx')\
                         .size().values
            xlabel = '%s Repeat Copies (exclusive)' % motif
            f = countPlot2(counts,target,xlabel,binsize=1)
            f.savefig(outfileName('.'.join([target,'motifCount_dist']),'png'),
                      dpi=DPI)
            #plot insert length
            counts = [insertLength[r] for r in data.readName.unique()]
            xlabel = 'Target Insert Length (bp)'
            f = countPlot2(counts,target,xlabel,binsize=1)
            f.savefig(outfileName('.'.join([target,'insertSize_dist']),'png'),
                      dpi=DPI)
        except KeyError:
            print '\tNo data for %s (skipping)' % target

    print 'Writing summary'
    summaryDf.to_csv(outfileName('summary','csv'))

    print 'Writing counts table'
    table = repeatDf.reset_index()\
                    .groupby(['target','readName','motif'])\
                    .size()\
                    .rename('repeats')\
                    .reset_index()
    writer = pd.ExcelWriter(outfileName('repeatCounts','xlsx'))
    for target,data in table.groupby('target'):
        motifs = [('repeats',m) for m in bed.query('name=="%s"'%target).motifs.iloc[0].split(',')]
        out = data[['readName','motif','repeats']].set_index(['motif','readName'])\
                                                  .unstack(0)\
                                                  .fillna(0).astype(int)\
                                                  .reindex(columns=motifs)\
                                                  .dropna(axis=1)
        #add in total repeat length
        out[(SIZECOL,'')] = out.index.map(insertLength)
        out.sort_values([SIZECOL,motifs[0]]).to_excel(writer,target)
    writer.close()

    return repeatDf,summaryDf

def getFlanks(ref,ctg,start,stop,flanksize=100,Lflank=None,Rflank=None):
    '''Extract flanking sequence from BED regions, given a pysam.FastaFile with .fai'''
    Lsize    = Lflank if Lflank else flanksize
    Rsize    = Rflank if Rflank else flanksize
    sequence = ref.fetch(ctg,start-Lsize,stop+Rsize)
    return [sequence[:Lsize],sequence[-Rsize:]]

def getFlankAligner(ref,ctg,start,stop,**kwargs):
    tmpRef = NamedTemporaryFile(mode='wb',delete=False)
    for side,seq in zip(['L','R'],getFlanks(ref,ctg,start,stop,**kwargs)):
        tmpRef.write('>{n}\n{s}\n'.format(n='_'.join([str(ctg),side]),s=seq))
    tmpRef.close()
    aligner = mp.Aligner(tmpRef.name)
    return aligner,tmpRef

class RepeatAnalysisReport_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='RepeatAnalysisReport.py', description='generate repeat-kde plot, waterfall plot, and summary for repeat expansions')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='input BAM of CCS alignments')
    parser.add_argument('inBED', metavar='inBED', type=str,
                    help='''input BED defining location of repeat region.
                            Repeats ONLY, no flanking sequence.\
                            Columns: ctg, start, end, name, motifs''')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta used for mapping BAM.  Must have .fai index.')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='Output directory.  Default %s' % os.getcwd())
    parser.add_argument('-s,--sample', dest='sample', type=str, default=None,
                    help='Sample name, prepended to output files.  Default None')
    parser.add_argument('-f,--flanksize', dest='flanksize', type=int, default=100,
                    help='Size of flanking sequence mapped for extracting repeat region.  Default 100')
    
    try:
        result = main(parser)
    except RepeatAnalysisReport_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
