#! /usr/bin/python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from pbcore.io import FastaReader,FastqReader

YLABEL     = 'CCS Reads'
DEFAULTBIN = 1
DEFAULTBW  = 3

def main(parser):
    args = parser.parse_args()
    
    fx            = args.inFastx if args.inFastx else sys.stdin
    #sortedRecs    = sorted(pysam.FastxFile(fx),key=lambda r:-len(r.sequence))
    try:
        recs = list(FastqReader(fx))
    except ValueError:
        #this will fail if fasta is streamed 
        recs = list(FastaReader(fx))

    counts = [rec.sequence.count(args.motif) for rec in recs]
    xlabel = f'{args.motif} Repeat Copies (exclusive)'
    labelMotif = list(map(eval,args.labelMotif.split(','))) if args.labelMotif else None
    f,c,b  = countPlot(counts,
                       args.name,
                       xlabel,
                       args.ylabel,
                       labelValues=labelMotif,
                       binsize=args.binsize,
                       bandwidth=args.bandwidth)
    f.savefig(f'{args.out}.motifcount.{args.format}',
              format=args.format,
              dpi=args.dpi)
    counts = list(map(len,recs))
    xlabel = 'Target Insert Length (bp)'
    labelLength = list(map(eval,args.labelLength.split(','))) if args.labelLength else None
    f,c,b  = countPlot(counts,
                       args.name,
                       xlabel,
                       args.ylabel,
                       labelValues=labelLength,
                       binsize=len(args.motif)*args.binsize,
                       bandwidth=len(args.motif)*args.bandwidth)
    f.savefig(f'{args.out}.insertSize.{args.format}',
              format=args.format,
              dpi=args.dpi)
    if args.exportBincounts:
        oname = f'{prgs.out}.histogramBins.csv'
        with open(oname,'w') as ofile:
            ofile.write('Length,Reads\n')
            for bn,cnt in zip(b,c):
                ofile.write(f'{bn},{cnt}\n')

    print('Done')
    return f

def countPlot(counts,title,xlabel,ylabel,labelValues=None,
              binsize=1,bandwidth=1):
    f,ax = plt.subplots()
    ax2 = ax.twinx()
    sns.kdeplot(counts,ax=ax2,color='k',bw=bandwidth,alpha=0.25)
    bincnts,bins,patches = ax.hist(counts,
                                   bins=range(min(counts)-1,
                                               max(counts)+1,
                                               binsize),
                                   align='left',
                                   alpha=0.5)
    maxY = max(bincnts)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.get_yaxis().set_visible(False)
    if labelValues:
        for val in labelValues:
            ax.annotate(str(val),
                        xy=(val,maxY), xycoords='data',
                        xytext=(val,1.1*maxY), textcoords='data',
                        arrowprops=dict(arrowstyle="->", color="0.25"),
                        horizontalalignment="center")
    return f,bincnts,bins

class CountPlot_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='plotCounts.py', description='generate histograms of motif counts and expansion size')
    parser.add_argument('-i','--inFastx', dest='inFastx', type=str, default=None, 
                    help='Input Fastx file. Default stdin')
    parser.add_argument('-o','--out', dest='out', type=str, default='hist',
                    help='Output prefix. default \'hist\'')
    parser.add_argument('-m','--motif', dest='motif', type=str, default=None, required=True,
                    help='Search motif')
    parser.add_argument('-b','--binsize', dest='binsize', type=int, default=DEFAULTBIN,
                    help='binsize for histogram. default %i' % DEFAULTBIN)
    parser.add_argument('-w','--bandwidth', dest='bandwidth', type=float, default=DEFAULTBW,
                    help='bandwidth for kde plot. default %f' % DEFAULTBW)
    parser.add_argument('-n','--name', dest='name', type=str, default='', required=False,
                    help='Title/name for figure')
    parser.add_argument('-y','--ylabel', dest='ylabel', type=str, default=YLABEL, required=False,
                    help='Y-axis label for figure.  Default %s' %YLABEL)
    parser.add_argument('-l','--labelMotif', dest='labelMotif', type=str, default=None,
                    help='Comma sep per-allele values to label with an arrow. Format \'2,35\'.  Default None')
    parser.add_argument('-L','--labelLength', dest='labelLength', type=str, default=None,
                    help='Comma sep per-allele values to label with an arrow. Format \'2,35\'.  Default None')
    parser.add_argument('-f','--format', dest='format', type=str, default='png',
                    help='Image format.  Default png')
    parser.add_argument('-d','--dpi', dest='dpi', type=int, default=400,
                    help='Image resolution.  Default 400')
    parser.add_argument('-e','--exportBins', dest='exportBincounts', action='store_true', default=False,
                    help='Export CSV of bin counts.  Default False')

    try:
        main(parser)
    except CountPlot_Exception as e:
        print('ERROR: %s' % e)
        sys.exit(1) 

