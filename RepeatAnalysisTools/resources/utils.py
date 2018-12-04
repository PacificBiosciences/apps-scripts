import pandas as pd
import numpy as np
import re

_RC_MAP = {'-': '-',
           'A': 'T',
           'C': 'G',
           'G': 'C',
           'N': 'N',
           'T': 'A',
           'a': 't',
           'c': 'g',
           'g': 'c',
           't': 'a'}

def rc(seq):
    '''revcomp'''
    return "".join([_RC_MAP[c] for c in seq[::-1]])

def getZmw(readname):
    '''extract the first two fields <movie>/<holenumber>/<other>'''
    return '/'.join(readname.split('/')[:2])

def readBED(bedfile,names=['ctg','start','end','name'],usecols=None):
    if not usecols:
        usecols = range(len(names))
    df = pd.read_table(bedfile,
                       sep='\s+',
                       names=names,
                       usecols=None)
    missing = df.columns[df.isna().apply(np.all) == True]
    if len(missing):
        raise ValueError('Missing columns in BED file: %s' % ','.join(map(str,missing)))
    return df

def getSubSeq(seq,aln):
    pos = sorted([getattr(a,att) for a in aln for att in ['q_st','q_en']])[1:-1]
    return seq[slice(*pos)]

def extractRepeat(sequence,aligner):
    aln = list(aligner.map(sequence))
    naln = len(aln)
    if naln == 2:
        seq = getSubSeq(sequence,aln)
        if aln[0].strand == -1:
            seq = rc(seq)
    else:
        seq = 'One Sided' if naln==1 else 'Poor/no Alignment'
    return seq

def getPositions(motif):
    def finder(seq):
        df = pd.DataFrame([1.0*m.start()/len(motif) for m in re.finditer(motif,seq)],
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
