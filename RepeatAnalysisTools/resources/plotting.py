import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import OrderedDict
import matplotlib.patches as mpatches

def waterfallPlot(data,row=None,x='position',y='idx',hue='motif'):
    g = sns.FacetGrid(data=data,
                      hue=hue,row=row,
                      height=5,aspect=2,
                      sharex=False,sharey=False)
    g.map(plt.scatter,x,y,marker=',',s=1.5)
    g.add_legend()
    g.set_xlabels('Position')
    g.set_ylabels('CCS Read')
    return g

def waterfallPlot2(data,row=None,x='position',y='idx',hue='motif'):
    motifLengths = data.motif.str.len()
    g = sns.FacetGrid(data=data,
                      hue=hue,row=row,
                      aspect=2,height=3,
                      sharex=False,sharey=False)
    for n in xrange(motifLengths.min()):
        g.data = data[motifLengths>=(n+1)].eval('position=position+@n')
        #g.map(sns.scatterplot,x,y,markers=[','],size=1.5)
        g.map(plt.scatter,x,y,marker=',',s=1.5)
    g.add_legend()
    g.set_xlabels('Position')
    g.set_ylabels('CCS Read')
    return g

def waterfallPlotRaster(data,x='position',y='idx'):
    palette = sns.color_palette()
    motifs  = data.motif.value_counts().index
    colors = OrderedDict([(m,palette[i]) for i,m in enumerate(motifs)])
    #all white grid
    raster = np.ones((int(data[y].max())+1,
                      int(data[x].max())+data.motif.str.len().max(),
                      3))
    #fill in motifs
    for i,row in data.iterrows():
        for j in xrange(len(row.motif)):
            raster[int(row[x]),int(row[y])+j,:] = colors[row.motif]
    #plot figure
    f,ax = plt.subplots()
    ax.imshow(raster,origin='lower',aspect='auto')
    ax.set_xlabel('Position')
    ax.set_ylabel('CCS Read')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    patches = [ mpatches.Patch(color=color, label=motif ) for motif,color in colors.items()]
    ax.legend(handles=patches,bbox_to_anchor=(1.25, 0.6),loc='best',frameon=False)
    
    return f

def countPlot(data,targetDict):
    '''targetdict is a dict of "target":"motif" listing primary motif to count'''
    name   = 'Repeat Copies'
    if 'target' not in data.columns:
        data['target'] = targetDict.keys()[0]
    #filter for primary motif
    qry    = '|'.join(['(target=="{t}" & motif=="{m}")'.format(t=t,m=ms.split(',')[0])
                       for t,ms in targetDict.items()])
    counts = data.query(qry)\
                 .groupby(['target','idx'])\
                 .size()\
                 .rename(name)\
                 .reset_index()
    g   = sns.FacetGrid(counts,
                        row='target',
                        height=4,
                        aspect=2,
                        sharex=False,
                        sharey=False)
    g.map(sns.kdeplot,name,bw=2)
    return g

def countPlot2(counts,target,xlabel,binsize=1):
    f,ax = plt.subplots()
    ax2 = ax.twinx()
    sns.kdeplot(counts,ax=ax2,color='k',bw=1,alpha=0.25)
    ax.hist(counts,
            bins=xrange(min(counts)-1,max(counts)+1,binsize),
            align='left',
            alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('CCS Reads')
    ax.set_title(target)
    ax.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.get_yaxis().set_visible(False)
    return f
