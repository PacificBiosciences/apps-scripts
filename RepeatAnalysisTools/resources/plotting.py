import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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

def countPlot2(counts,target,motif,binsize=1):
    f,ax = plt.subplots()
    ax.hist(counts,
            bins=xrange(min(counts)-1,max(counts)+1,binsize),
            align='left',
            alpha=0.6)
    ax.set_xlabel(' '.join([motif,'Repeat Copies']))
    ax.set_ylabel('CCS Reads')
    ax.set_title(target)
    ax2 = ax.twinx()
    sns.kdeplot(counts,ax=ax2,color='k',bw=1,alpha=0.25)
    ax.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.get_yaxis().set_visible(False)
    return f
