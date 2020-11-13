#! /usr/bin/python

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
import sys,re,pysam
from operator import attrgetter

COLORMAP = cm.get_cmap('tab10')
UNKNOWN  = (0.75,)*3
BLANK    = (1,)*3
XLABEL   = 'Position'
YLABEL   = 'CCS Read'
MAXQV    = 93
MINQV    = 0
CENTERQV = 20
QVCOLOR  = cm.get_cmap('RdYlBu')

def main(parser):
    args = parser.parse_args()
    
    infastx    = args.inFastx if args.inFastx else "-"
    motifs     = args.motifs.split(',')
    sortedRecs = sorted(list(pysam.FastxFile(infastx,'r')),key=sortFunc(args.sortCluster))

    colors  = OrderedDict([(m,COLORMAP.colors[i]) for i,m in enumerate(motifs)])
    raster  = motifRaster(sortedRecs,motifs,colors)
    patches = [ mpatches.Patch(color=color, label=motif ) for motif,color in list(colors.items())]
    f,ax    = plotWaterfall(raster,XLABEL,args.ylabel,labels=patches)

    out = args.out if args.out.endswith(args.format) else '%s.%s' % (args.out,args.format) 
    plt.tight_layout()
    f.savefig(out,dpi=args.dpi,format=args.format)

    if args.plotQV:
        qvraster = qvRaster(sortedRecs,'Base QV')
        norm     = TwoSlopeNorm(CENTERQV,vmin=MINQV,vmax=MAXQV)
        f,ax     = plotWaterfall(qvraster,XLABEL,args.ylabel,norm=norm,cmap=QVCOLOR,colorbar='QV')
        name,ext = out.rsplit('.',1)
        f.savefig(f'{name}.QV.{ext}',format=args.format)
    
    print('Done')
    return raster

def sortFunc(csvMap):
    if csvMap:
        clusterMap = dict(s.split(',') for s in open(csvMap).read().split())
        def getVal(rec):
            name = '/'.join(rec.name.split('/')[:3])
            try:
                return int(clusterMap.get(name,-1)),-len(rec.sequence)
            except ValueError:
                raise Waterfall_Exception('Clusters must be numeric') 
        return getVal
    else:
        return lambda rec: -len(rec.sequence)

def plotWaterfall(array,xlabel,ylabel,labels=None,colorbar=False,**kwargs):
    f,ax  = plt.subplots()
    image = ax.imshow(array,origin='lower',aspect='auto',**kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if labels:
        ax.legend(handles=labels,bbox_to_anchor=(1.25, 0.6),loc='best',frameon=False)
    if colorbar:
        f.colorbar(image,shrink=0.6, label=colorbar)
    return f,ax
    
def motifRaster(recs,motifs,colors):
    raster        = np.ones((len(recs),
                             len(recs[0].sequence),
                             3))
    patt          = re.compile('|'.join(['(%s)'%m for m in motifs]))
    colors['other'] = UNKNOWN

    for i,rec in enumerate(recs):
        for j in patt.finditer(rec.sequence):
            raster[i,j.start():j.end(),:] = colors[j.group()]

        #fill in unknown cells
        blank = np.all(raster[i] == BLANK,axis=1)
        blank[len(rec.sequence):] = False
        raster[i,blank,:] = colors['other']
    return raster

def qvRaster(recs,title):
    raster  = -np.ones((len(recs),len(recs[0].sequence)))
    for i,rec in enumerate(recs):
        raster[i,:len(rec.sequence)] = phred2QV(rec)
    return np.ma.masked_where(raster == -1,raster)

def phred2QV(rec):
    return np.array([ord(q)-33 for q in rec.quality])

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
    parser.add_argument('-s,--sortCluster', dest='sortCluster', type=str, default=None,
                    help='Sort reads by cluster, determined by file of <readname>,<cluster>. Default sort by length only')
    parser.add_argument('-y,--ylabel', dest='ylabel', type=str, default=YLABEL, required=False,
                    help='Y-axis label for plot. Default %s'%YLABEL)
    parser.add_argument('-f,--format', dest='format', type=str, default='png',
                    help='Image format.  Default png')
    parser.add_argument('-d,--dpi', dest='dpi', type=int, default=400,
                    help='Image resolution.  Default 400')
    parser.add_argument('-q,--plotQV', dest='plotQV', action='store_true', default=False,
                    help='Plot additional QV waterfall.  Default False')

    try:
        main(parser)
    except Waterfall_Exception as e:
        print(f'ERROR:{e}')
        sys.exit(1) 

