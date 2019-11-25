#! /usr/bin/python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import matplotlib.patches as mpatches
from matplotlib import cm
import sys,re
from pbcore.io import FastaReader,FastqReader

COLORMAP = cm.get_cmap('tab10')
UNKNOWN  = (0.75,)*3
BLANK    = (1,)*3
XLABEL   = 'Position'
YLABEL   = 'CCS Read'

def main(parser):
    args = parser.parse_args()
    
    fx            = args.inFastx if args.inFastx else sys.stdin
    #sortedRecs    = sorted(pysam.FastxFile(fx),key=lambda r:-len(r.sequence))
    try:
        sortedRecs = sorted(FastqReader(fx),key=len,reverse=True)
    except ValueError:
        #this will fail if fasta is streamed 
        sortedRecs = sorted(FastaReader(fx),key=len,reverse=True)

    raster        = np.ones((len(sortedRecs),
                             len(sortedRecs[0].sequence),
                             3))
    motifs        = args.motifs.split(',')
    colors        = OrderedDict([(m,COLORMAP.colors[i]) for i,m in enumerate(motifs)])
    colors['other'] = UNKNOWN

    for i,rec in enumerate(sortedRecs):
        ##iter through backwards so any 
        ##greedy match over-writes end up 
        ##as the first listed motif 
        #for motif in motifs[::-1]:
        #    for j in re.finditer(motif,rec.sequence):
        #        raster[i,j.start():j.end(),:] = colors[motif]
        #better, search for all at same time
        patt = re.compile('(%s)' % (')|('.join(motifs)))
        for j in patt.finditer(rec.sequence):
            raster[i,j.start():j.end(),:] = colors[j.group()]

        #fill in unknown cells
        blank = np.all(raster[i] == BLANK,axis=1)
        blank[len(rec.sequence):] = False
        raster[i,blank,:] = colors['other']

    f,ax = plt.subplots()
    ax.imshow(raster,origin='lower',aspect='auto')
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    patches = [ mpatches.Patch(color=color, label=motif ) for motif,color in colors.items()]
    ax.legend(handles=patches,bbox_to_anchor=(1.25, 0.6),loc='best',frameon=False)

    out = args.out if args.out.endswith(args.format) else '%s.%s' % (args.out,args.format) 
    plt.tight_layout()
    f.savefig(out,dpi=args.dpi,format=args.format)
    print 'Done'
    return f

class Waterfall_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='waterfall.py', description='quick waterfall plot from fastx of extracted repeat sequences')
    parser.add_argument('-i,--inFastx', dest='inFastx', type=str, default=None, 
                    help='Input Fastx file. Default stdin')
    parser.add_argument('-o,--out', dest='out', type=str, default=None, required=True,
                    help='Output file')
    parser.add_argument('-m,--motifs', dest='motifs', type=str, default=None, required=True,
                    help='Search motifs, comma separated, most frequent first, e.g. \'CGG,AGG\'')
    parser.add_argument('-f,--format', dest='format', type=str, default='png',
                    help='Image format.  Default png')
    parser.add_argument('-d,--dpi', dest='dpi', type=int, default=400,
                    help='Image resolution.  Default 400')

    try:
        main(parser)
    except Waterfall_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1) 

