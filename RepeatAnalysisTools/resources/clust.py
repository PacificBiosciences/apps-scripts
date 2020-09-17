import pysam,re
import numpy as np
import pandas as pd
from collections import Counter
from operator import itemgetter

NOCLUST  = 999
DEFAULTCI= [0.025,0.975] #95%

def addHPtag(inBAM,outBAM,clusterMap,noCluster=NOCLUST,dropNoClust=False):
    '''clusterMap is map of {readname:cluster::int}'''
    with pysam.AlignmentFile(inBAM) as inbam:
        with pysam.AlignmentFile(outBAM,'wb',template=inbam) as outbam:
            for rec in inbam:
                if rec.query_name not in clusterMap and dropNoClust:
                    continue 
                clust = int(clusterMap.get(rec.query_name,noCluster))
                rec.set_tag('HP',clust)
                outbam.write(rec)
    pysam.index(outBAM)
    return None

def clusterName(vals):
    'vals: (clusterInt,numreads)'
    return 'cluster{0[0]}_numreads{0[1]}'.format(vals)

def getCluster(name):
    return int(re.search('cluster(\d+)',name).group(1))

def readClusterFile(clustfile,nFields=3):
    '''
    nFields: use n /-separated fields from input read names in result.
             nFields = 0 -> all fields
    '''
    res  = {}
    name,cluster = None,None
    with open(clustfile) as f:
        for line in f:
            if line.startswith('>'):
                cluster = getCluster(line[1:])
            else:
                read = '/'.join(line.split('/')[:nFields]) if nFields else line.strip()
                res[read] = cluster
    return res

def getCounts(seqGen,kmerSize,motifCounter):
    counts = {name:[pd.Series(getKmerCounts(seq,k=kmerSize)),
                    pd.Series(motifCounter(seq))]
              for name,seq,qual in seqGen}
    kmer,motif = [pd.DataFrame(list(map(itemgetter(i),list(counts.values()))),
                               index=list(counts.keys())).fillna(0)
                  for i in [0,1]]
    return kmer,motif

def getKmerCounts(seq,k=3):
    return Counter(seq[i:i+k] for i in range(0,len(seq)-k+1))

def resampleCI(data,nboot=10000,ci=DEFAULTCI):
    n = max(nboot,len(data))
    resamp = np.random.choice(data,size=n,replace=True)
    return '({} - {})'.format(*list(map(int,np.quantile(resamp,ci))))

def clusterStats(motifCounts,clusterIdx,outColumns,
                 aggFuncs=[np.median,np.mean],
                 randomSeed=None,ci=DEFAULTCI):
    '''
    motifCounts: df with cols =  motif (+length), index = readnames
    clusterIdx: vector of cluster indices, same order as motifCounts.index
    outColumns: list of column names to describe in output
    aggFuncs: list of functions to apply to each column
    randomSeed: random seed for resampling ci
    '''
    clusters    = motifCounts.groupby(clusterIdx)
    clusterSize = clusters.size().rename(('Read','count'))

    #set random seed
    np.random.seed(randomSeed)
    results = clusters[outColumns].agg(aggFuncs+[resampleCI])\
                                  .join(clusterSize)
    #rename clusters
    names = clusterSize.reset_index().apply(clusterName,axis=1)
    results.index = names.values

    #rename column
    name = 'ci%i' % int(100*(ci[1]-ci[0]))
    results.rename(columns={'resampleCI':name},level=1,inplace=True)    

    return names,results
