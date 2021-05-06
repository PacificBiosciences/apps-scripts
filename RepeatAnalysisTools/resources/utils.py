import pandas as pd
import numpy as np
import re
from collections import Counter,OrderedDict
from operator import itemgetter
from itertools import groupby

def readBED(bedfile,names=['ctg','start','end','name'],usecols=None):
    if not usecols:
        usecols = list(range(len(names)))
    dtypes = {'ctg':np.str,'start':np.int64,'stop':np.int64}
    df = pd.read_csv(bedfile,
                     sep='\s+',
                     names=names,
                     usecols=usecols,
                     dtype=dtypes)
    missing = df.columns[df.isna().all()]
    if len(missing):
        raise RepeatAnalysisUtils_Exception('Missing columns in BED file')#: %s' % ','.join(map(str,missing)))
    return df

def countMotifs(motifs,lengthField=False,collapseHP=False,blocks=False):
    ''' takes a list of motifs
        returns a function that takes a sequence and returns a dict of counts
        optionally includes length of seq in output
        optionally hpcollapse motifs and reads before counting
        **length is always original sequence, not hpcollapsed**
        if "blocks" counts blocks of first motif sep by second
    '''
    transform =  hpCollapse if collapseHP else (lambda x:x)
    mm       = OrderedDict(list(zip(list(map(transform,motifs)),motifs)))
    if len(mm) < len(motifs): 
        vals = motifs.__repr__(),list(map(hpCollapse,motifs)).__repr__()
        e    = 'Degenerate collapsed motifs\noriginal {}\ncollapsed {}\n'.format(*vals)
        raise RepeatAnalysisUtils_Exception(e)
    patt      = re.compile('(' + ')|('.join(list(mm.keys())) + ')')
    blockpatt = re.compile('(' + ')|('.join(list(mm.keys())[:2]) + ')') 
    def getCounts(seq):
        counts = Counter(mm[m.group()] for m in patt.finditer(transform(seq)))
        #fill for motifs not found 
        counts.update({m:0 for m in list(mm.values()) if m not in counts})
        if lengthField:
            counts[lengthField] = len(seq)
        if blocks:
            blocked = tuple((m,len(list(g))) for m,g in groupby([mm[m.group()] 
                       for m in blockpatt.finditer(transform(seq))]))
            counts[blocks] = counts2name(blocked,motifs[0])
        return counts
    return getCounts

def counts2name(cntList,primary=None):
    tot  = sum(map(itemgetter(1),cntList))
    name = '+'.join(str(c) for m,c in cntList if m==primary or primary is None)
    return f'{tot}:{name}'

_PATT = re.compile(r'([ATGC])\1+')
def hpCollapse(seq):
    return _PATT.sub(r'\1',seq)

def formatMixFloat(fmt):
    def formatter(elem):
        if pd.isnull(elem):
            return ''
        try: #integers not in float fmt
            if int(elem) == elem:
                return elem
            else:
                return fmt % elem
        except (TypeError,ValueError):
            return elem
    return formatter

class RepeatAnalysisUtils_Exception(Exception):
    pass
