import pandas as pd
import numpy as np
import mappy as mp
import re
from textwrap import wrap
from tempfile import NamedTemporaryFile

_RC_MAP = dict(zip('-ACGNTacgt','-TGCNAtgca'))

def rc(seq):
    '''revcomp'''
    return "".join([_RC_MAP[c] for c in seq[::-1]])

def getZmw(readname):
    '''extract the first two fields <movie>/<holenumber>/<other>'''
    return '/'.join(readname.split('/')[:2])

def getHn(name):
    '''extract the holenumber from a readname'''
    return int(name.split('/')[1])

def readBED(bedfile,names=['ctg','start','end','name'],usecols=None):
    if not usecols:
        usecols = range(len(names))
    dtypes = {'ctg':np.str,'start':np.int64,'stop':np.int64}
    df = pd.read_table(bedfile,
                       sep='\s+',
                       names=names,
                       usecols=usecols,
                       dtype=dtypes)
    missing = df.columns[df.isna().apply(np.all) == True]
    if len(missing):
        raise RepeatAnalysisUtils_Exception('Missing columns in BED file')#: %s' % ','.join(map(str,missing)))
    return df

def getSubSeq(seq,aln):
    pos = sorted([getattr(a,att) for a in aln for att in ['q_st','q_en']])[1:-1]
    return pos + [seq[slice(*pos)]]
    #return seq[slice(*pos)]

def extractRepeat(sequence,aligner):
    aln = list(aligner.map(sequence))
    naln = len(aln)
    if naln == 2:
        start,stop,seq = getSubSeq(sequence,aln)
        #seq = getSubSeq(sequence,aln)
        if aln[0].strand == -1:
            seq = rc(seq)
    else:
        seq = 'One Sided' if naln==1 else 'Poor/no Alignment'
        start,stop = None,None
    return start,stop,seq
    #return seq

def getPositions(motif):
    def finder(seq):
        #df = pd.DataFrame([1.0*m.start()/len(motif) for m in re.finditer(motif,seq)],
        #                   columns=['position'])
        df = pd.DataFrame([1.0*m.start() for m in re.finditer(motif,seq)],
                           columns=['position'])
        df['motif'] = motif
        return df
    return finder

def countAlignments(repeatRegions,reference=None):
    summary = pd.Series({'totalReads': len(repeatRegions)})
    oneSided = repeatRegions.query('subsequence=="One Sided"').index
    poorAln  = repeatRegions.query('subsequence=="Poor/no Alignment"').index
    summary['oneSided']      = len(oneSided)
    summary['poorAlignment'] = len(poorAln)
    summary['spanningReads'] = len(repeatRegions) - (summary['oneSided'] + summary['poorAlignment'])
    if reference:
        summary['reference'] = reference
    return summary,repeatRegions.drop(oneSided.append(poorAln))

def writeFasta(fastaName,records,nWrap=60):
    with open(fastaName,'w') as fa:
        for name,sequence in records:
            seq = '\n'.join(wrap(sequence,nWrap))
            fa.write('>{name}\n{seq}\n'.format(name=name,seq=seq))
    return fastaName

def getFlanks(ref,ctg,start,stop,flanksize=100,Lflank=None,Rflank=None):
    '''Extract flanking sequence from BED regions, given a pysam.FastaFile with .fai'''
    Lsize    = Lflank if Lflank else flanksize
    Rsize    = Rflank if Rflank else flanksize
    sequence = ref.fetch(ctg,start-Lsize,stop+Rsize)
    return [sequence[:Lsize],sequence[-Rsize:]]

def getFlankAligner(ref,ctg,start,stop,**kwargs):
    tmpRef = NamedTemporaryFile(mode='wb',delete=False)
    for side,seq in zip(['L','R'],getFlanks(ref,ctg,start,stop,**kwargs)):
        tmpRef.write('>{n}\n{s}\n'.format(n='_'.join([str(ctg),side]),s=seq))
    tmpRef.close()
    aligner = mp.Aligner(tmpRef.name)
    return aligner,tmpRef

class RepeatAnalysisUtils_Exception(Exception):
    pass
