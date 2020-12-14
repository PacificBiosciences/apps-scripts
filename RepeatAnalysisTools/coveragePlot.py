#! /usr/bin/python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pysam,os
import pandas as pd
import numpy as np
import seaborn as sns
from resources.utils import readBED

DWINDOW  = 2500 
FIGHEIGHT= 4
FIGWIDTH = 16
XLABEL   = 'Chromosome'
YLABEL   = 'Reads'
PALETTE  = 'husl'
DPI      = 400

#Define some human chromosome names
#add more sets as needed
#all and only these names will plot
withCH = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
          'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',
          'chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM']
CHROM    = {'hs37d5': ['1','2','3','4','5','6','7','8','9','10','11','12','13',
                       '14','15','16','17','18','19','20','21','22','X','Y','MT'],
            'hg19'  : withCH,
            'hg38'  : withCH 
           }            
def main(parser):
    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.inBAM)

    if args.nameFormat == 'all':
        useChroms = bam.references
    elif args.nameFormat is None:
        useChroms = None
        print(f'Guessing chromosome name format')
        for name,chroms in list(CHROM.items()):
            if chroms[0] in bam.references:
                useChroms = chroms
                break
    else:
        useChroms = CHROM[args.nameFormat]
    if not useChroms:
        raise CoveragePlot_Exception('Chromosome set not found. Check CHROM definition in script')
    print(f'Using {args.nameFormat if args.nameFormat else name} chromosome name format')

    #get all covered positions in a df
    #recommended only for sparse (targeted) coverage!
    print('Reading coverage')
    cov = pd.DataFrame([(chrom,pile.pos,pile.nsegments)
                        for chrom in bam.references
                        for pile in bam.pileup(chrom,stepper='all')
                        if chrom in useChroms],
                       columns=['chr','pos','coverage'])

    #map to get nbins given window size and contig length
    nbins     = {chrom : length//args.window + 1
                 for chrom,length in zip(bam.references,bam.lengths)}
    #vector of intervals
    intervals = pd.cut(cov.pos,range(0,max(bam.lengths),args.window))
    #average for each window
    print('Calculating mean values')
    try:
        meancov   = cov.groupby(['chr',intervals]).coverage.mean()
    except:
        raise CoveragePlot_Exception(f'No coverage found for any of contigs {useChroms}')
    #index of intervals for ordering results
    cats      = intervals.cat.categories

    #if targets BED, set index to interval bin number
    if args.targets:
        targets = readBED(args.targets)
        targets.index = pd.cut(targets.eval('(start+end)/2'),cats).cat.codes
    #else set targets to empty df
    else:
        targets = pd.DataFrame(columns=['ctg'])

    #plot
    fig,ax = plt.subplots()
    fig.set_figheight(FIGHEIGHT)
    fig.set_figwidth(FIGWIDTH)
    maxy   = meancov.max()
    xstart = 0 
    label  = 1
    ticks,ticklabels = [],[]
    with sns.color_palette(PALETTE, len(useChroms)):
        for chrom in useChroms:
            try:
                bins = nbins[chrom]
                print(f'Plotting coverage: {chrom}')
            except KeyError as e:
                print(f'Contig {chrom} not found, skipping')
                continue
            if chrom in meancov:
                #order values and fill where no coverage
                yvals = meancov[chrom].reindex(cats[:bins],fill_value=0).fillna(0)
            else: #no coverage at all
                yvals = np.zeros(bins)
            xstop = xstart+len(yvals)
            xvals = range(xstart,xstop)
            p     = ax.plot(xvals,yvals)[0]
            if label: #every other chrom
                ticks.append(xstart + len(yvals)/2)
                ticklabels.append(chrom)
            #label any targets in the chrom
            for pos,target in targets.query('ctg==@chrom').iterrows():
                xpos = xstart+pos
                annotateTarget(ax,target['name'],xpos,maxy,p.get_color())
            label^=1  #flip
            xstart = xstop
            ax.vlines(xstart,0,maxy,colors='k',alpha=0.25,linestyles=':')
    #esthetics
    for loc in ['top','bottom','right']:
        ax.spines[loc].set_visible(False)
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    ax.xaxis.set_tick_params(length=0)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)

    #save figure
    name = f'{args.prefix}{"/" if os.path.isdir(args.prefix) else "."}coveragePlot.png'
    fig.savefig(name,dpi=DPI)
    print(f'Saved {name}')

def annotateTarget(ax,text,xpos,maxy,color):
    ax.annotate(text, xy=(xpos,1.01*maxy), xytext=(xpos,1.11*maxy),
                arrowprops=dict(facecolor=color,edgecolor=color,
                                width=1,headwidth=7,
                                headlength=9,shrink=0.001),
                horizontalalignment='center')
    return None

class CoveragePlot_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse,sys

    parser = argparse.ArgumentParser(prog='coveragePlot.py', description='Plot coverage across reference')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='BAM file of aligned reads.  Must have .bai index')
    parser.add_argument('-o,--outprefix', dest='prefix', type=str, default=os.getcwd(),
                    help='prefix for output file.  default cwd.')
    parser.add_argument('-w,--window', dest='window', type=int, default=DWINDOW,
                    help=f'window size for averaging coverage across reference.  default {DWINDOW}')
    parser.add_argument('-t,--targetBED', dest='targets', type=str, default=None,
                    help='label plot with targets from BED file.  default None')
    parser.add_argument('-n,--nameFormat', dest='nameFormat', choices=list(CHROM.keys()) + ['all'], default=None,
                    help='''Chromosome name format. 
                          Limits plotted set to 23 chrom and mito.
                          hg19/hg38: "chr1", "chr2",...
                          hs37d5: 1,2,...
                          Use "all" to use all mapped references. default guess''')

    try:
        main(parser)
    except CoveragePlot_Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1) 





