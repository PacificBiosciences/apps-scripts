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
    g.set_xlabels('Repeat Copies')
    g.set_ylabels('CCS Read')
    return g

def countPlot(data,targetDict):
    '''targetdict is a dict of "target":"motif" listing primary motif to count'''
    name   = 'Repeat Copies'
    if 'target' not in data.columns:
        data['target'] = targetDict.keys()[0]
    #filter for primary motif
    qry    = '|'.join(['(target=="%s" & motif=="%s")' % tup for tup in targetDict.items()])
    counts = data.query(qry)\
                 .groupby(['target','idx'])\
                 .size()\
                 .rename(name)\
                 .reset_index()
    g   = sns.FacetGrid(counts,row='target',height=4,aspect=2)
    g.map(sns.kdeplot,name,bw=2)
    return g

