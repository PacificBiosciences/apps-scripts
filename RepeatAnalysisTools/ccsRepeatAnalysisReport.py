#! /home/UNIXHOME/jharting/anaconda2/bin/python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import re,os,sys
from pbcore.io import FastaReader,FastqReader
import mappy as mp

sns.set_style('white')
DPI=400
FMR1FLANKS='resources/FMR1_L446_R503.fasta'

_RC_MAP = {'-': '-',
           'A': 'T',
           'C': 'G',
           'G': 'C',
           'N': 'N',
           'T': 'A',
           'a': 't',
           'c': 'g',
           'g': 'c',
           't': 'a'}

repeatPatterns = {'FMR1' : ['CGG','AGG'],
                  'HTT'  : ['CAG'],
                  'ALS'  : ['GGCCCC'],
                  'FUCHS': ['TGC']}

def main(parser):
    args = parser.parse_args()

    if args.motifs:
        print 'Over-riding preset motifs with %s' % args.motifs
        motifs = args.motifs.split(',')
        label  = args.label
    else:
        motifs = repeatPatterns[args.preset]
        label  = args.label if args.label else args.preset

    reader  = FastaReader if args.ccsFastx.endswith('a') else FastqReader
    aligner = mp.Aligner(args.target)
    summary = pd.Series(index=['totalReads',
                               'spanningReads',
                               'oneSided',
                               'poorAlignment',
                               'reference'])
    #                           'repeat_est1',
    #                           'repeat_est2'])
    summary['reference'] = args.target

    #function to generate output names
    sample      = args.sample + '.' if args.sample else ''
    label       = label + '.' if label else ''
    outfileName = lambda name,ext: '{d}/{s}{l}{n}.{e}'.format(d=args.outDir,s=sample,l=label,n=name,e=ext)

    #function to write summary
    writeSummary = lambda: summary.to_csv(outfileName('summary','csv'))

    print 'Mapping and extracting repeat regions'
    repeatRegions = pd.DataFrame({'read':rec.name, 
                                  'subsequence': extractRepeat(rec.sequence,aligner)}
                                 for rec in reader(args.ccsFastx))

    #count before and after removing reads that do not align or
    summary['totalReads']    = len(repeatRegions)
    oneSided = repeatRegions.query('subsequence=="One Sided"').index
    summary['oneSided']      = len(oneSided)
    poorAln  = repeatRegions.query('subsequence=="Poor/no Alignment"').index 
    summary['poorAlignment'] = len(poorAln)
    for labels in [oneSided,poorAln]:
        repeatRegions.drop(labels,inplace=True)
    summary['spanningReads'] = len(repeatRegions)

    repeatRegions = repeatRegions.assign(size=repeatRegions.subsequence.map(len))\
                                 .sort_values('size',ascending=False)\
                                 .drop(columns='size')\
                                 .reset_index(drop=True)    

    print 'Counting repeats'
    try:
        motifDfs = [pd.concat(repeatRegions.set_index('read',append=True).subsequence.map(getPositions(motif)).to_dict())\
                      .reset_index(level=2,drop=True)\
                      .reset_index()\
                      .rename(columns={'level_0':'idx','level_1':'readName'})
                    for motif in motifs]
    except ValueError,e:
        writeSummary()
        raise ccsWaterFallRepeats_Exception('No reads map to target!')

    allDf = pd.concat(motifDfs,ignore_index=True)

    print 'Plotting Waterfall'
    g = sns.FacetGrid(data=allDf,hue='motif',height=6,aspect=2)
    g.map(plt.scatter,'position','idx',marker=',',s=2)
    g.add_legend()
    g.set_xlabels('Repeat Copies')
    g.set_ylabels('CCS Read')
    g.savefig(outfileName('waterfall','png'),dpi=DPI)

    print 'Plotting histogram'
    name = 'Repeat Copies'
    data = allDf.query('motif=="%s"' % motifs[0])\
                .groupby('idx')\
                .size()\
                .rename(name)\
                .to_frame()
    g = sns.FacetGrid(data,height=4,aspect=2)
    g.map(sns.kdeplot,name,bw=1)
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

def rc(seq):
    return "".join([_RC_MAP[c] for c in seq[::-1]])

def getSubSeq(seq,aln):
    pos = sorted([getattr(a,att) for a in aln for att in ['q_st','q_en']])[1:-1]
    return seq[slice(*pos)]

def extractRepeat(sequence,aligner):
    aln = list(aligner.map(sequence))
    naln = len(aln)
    if naln == 2:
        seq = getSubSeq(sequence,aln)
        if aln[0].strand == -1:
            seq = rc(seq)
    else:
        seq = 'One Sided' if naln==1 else 'Poor/no Alignment'
    return seq

def getPositions(motif):
    def finder(seq):
        df = pd.DataFrame([m.start()/3. for m in re.finditer(motif,seq)],
                           columns=['position'])
        df['motif'] = motif
        return df
    return finder

class ccsWaterFallRepeats_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='ccsRepeatAnalysisReport.py', description='generate repeat-kde plot,waterfall plot, and summary for repeat expansion ccs reads')
    parser.add_argument('ccsFastx', metavar='ccsFastx', type=str,
                    help='input ccs reads, can be fasta or fastq')
    parser.add_argument('-t,--target', dest='target', type=str, 
                    default=FMR1FLANKS,
                    help='fasta reference of sequence flanking repeat region (2 sequences). Default %s' % FMR1FLANKS)
    parser.add_argument('-p,--preset', dest='preset', type=str, default='FMR1',
                    help='preset motifs to search for.  default \'FMR1\'.  Available: %s' % ','.join(repeatPatterns.keys()))
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
    except ccsWaterFallRepeats_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
