#! /home/UNIXHOME/jharting/anaconda2/bin/python

import numpy as np
import pandas as pd
import os,sys
import mappy as mp
from resources.utils import extractRepeat,\
                            countAlignments,\
                            getPositions
from resources.plotting import waterfallPlot,countPlot

DPI=400
FMR1FLANKS='%s/resources/FMR1_L446_R503.fasta' % sys.path[0]

repeatPatterns = {'FMR1' : ['CGG','AGG'],
                  'HTT'  : ['CAG','CAA', 'CCG', 'CCA', 'CGG'],
                  'ALS'  : ['GGCCCC'],
                  'FUCHS': ['TGC'],
                  'Sca10': ['ATTCT', 'ATTCC', 'ATTTCT', 'ATTCTT']}

def main(parser):
    args = parser.parse_args()

    if args.motifs:
        print 'Over-riding preset motifs with %s' % args.motifs
        motifs = args.motifs.split(',')
        label  = args.label
    else:
        motifs = repeatPatterns[args.preset]
        label  = args.label if args.label else args.preset

    aligner = mp.Aligner(args.target)

    #function to generate output names
    s           = args.sample + '.' if args.sample else ''
    l           = label + '.' if label else ''
    outfileName = lambda name,ext: '{d}/{s}{l}{n}.{e}'.format(d=args.outDir,s=s,l=l,n=name,e=ext)

    #function to write summary
    writeSummary = lambda: summary.to_csv(outfileName('summary','csv'))

    print 'Mapping and extracting repeat regions'
    repeatRegions = pd.DataFrame({'read':rec[0], 
                                  'subsequence': extractRepeat(rec[1],aligner)}
                                 for rec in mp.fastx_read(args.ccsFastx))

    repeatRegions = repeatRegions.assign(size=repeatRegions.subsequence.map(len))\
                                 .sort_values('size',ascending=False)\
                                 .drop(columns='size')\
                                 .reset_index(drop=True)    
    #filter and summarize
    summary,filtered = countAlignments(repeatRegions,reference=args.target)


    print 'Counting repeats'
    try:
        motifDfs = [pd.concat(filtered.set_index('read',append=True).subsequence.map(getPositions(motif)).to_dict())\
                      .reset_index(level=2,drop=True)\
                      .reset_index()\
                      .rename(columns={'level_0':'idx','level_1':'readName'})
                    for motif in motifs]
    except ValueError,e:
        writeSummary()
        raise fastRepeatAnalysisReport_Exception('No reads map to target!')

    allDf = pd.concat(motifDfs,ignore_index=True)

    print 'Plotting Waterfall'
    g = waterfallPlot(allDf)
    g.savefig(outfileName('waterfall','png'),dpi=DPI)

    print 'Plotting histogram'
    g = countPlot(allDf,{label:motifs[0]})
    g.savefig(outfileName('repeatCount_kde','png'),dpi=DPI)

    print 'Writing summary'
    writeSummary()

    print 'Writing counts table'
    allDf.groupby(['readName','motif'])\
         .size()\
         .unstack()\
         .fillna(0)\
         .astype(int)\
         .to_csv(outfileName('repeatCounts','csv'))

    return allDf

class fastRepeatAnalysisReport_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='fastxRepeatAnalysisReport.py', description='generate repeat-kde plot,waterfall plot, and summary for repeat expansion ccs reads')
    parser.add_argument('ccsFastx', metavar='ccsFastx', type=str,
                    help='input ccs reads, can be fasta or fastq')
    parser.add_argument('-t,--target', dest='target', type=str, 
                    default=FMR1FLANKS,
                    help='fasta reference of sequence flanking repeat region (2 sequences). Default %s' % FMR1FLANKS)
    parser.add_argument('-p,--preset', dest='preset', type=str, default='FMR1',
                    help='preset motifs to search for.  default \'FMR1\'.\nAvailable: %s' % ','.join(repeatPatterns.keys()))
    parser.add_argument('-m,--motifs', dest='motifs', type=str, default=None,
                    help='Comma separated motifs to search for.  Will over-ride presets.  Default none.')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='Output directory.  Default %s' % os.getcwd())
    parser.add_argument('-s,--sample', dest='sample', type=str, default=None,
                    help='Sample name, prepended to output files.  Default None')
    parser.add_argument('-l,--label', dest='label', type=str, default=None,
                    help='Label, prepended to output files after sample name, e.g. \'FMR1\'.  Defaults to preset name, except if -m set')
    
    try:
        result = main(parser)
    except fastRepeatAnalysisReport_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
