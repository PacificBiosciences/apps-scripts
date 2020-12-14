import os,re,pysam
import mappy as mp
import numpy as np
from tempfile import NamedTemporaryFile

ALIGNFILTER=0x900
MINOUTPUTLEN=6

def extractRegion(inBAM,reference,region=None,ctg=None,start=None,stop=None,flanksize=100,revcomp=False,minRQ=None):
    ref = pysam.FastaFile(reference)
    bam = pysam.AlignmentFile(inBAM)
    if region:
        try:
            ctg,start,stop = getCoordinates(region)
        except AttributeError:
            #catch when the region format doesn't match
            raise Extract_Exception(f'Invalid region format {region}. Correct \'[chr]:[start]-[stop]\'')
    elif ctg==None or start==None or stop==None:
            #catch when missing coord and no region passed
            raise Extract_Exception('Must pass either valid region string or all of ctg,start,stop')

    aligner,tmp = getFlankAligner(ref,ctg,start-1,stop,flanksize=flanksize) 
    
    try:
        for rec in bam.fetch(ctg,start,stop):
            if (rec.flag & ALIGNFILTER):
                continue
            rStart,rStop,subseq = extractRepeat(rec.query_sequence,aligner)
            if len(subseq) <= MINOUTPUTLEN:
                continue
            if rStart:
                name   = nameFunction(rec.query_name,rStart,rStop)
                qvVals = rec.query_qualities[rStart:rStop] 
                if minRQ:
                    meanrq = meanRQ(qvVals)
                    if meanrq < minRQ:
                        continue
                qual = ''.join([chr(q+33) for q in qvVals])
                if revcomp:
                    subseq = rc(subseq)
                    qual   = qual[::-1]
                yield name,subseq,qual
    finally:
        os.remove(tmp.name)

def fqRec(name,seq,qual):
    return f'@{name}\n{seq}\n+\n{qual}\n'

def getCoordinates(regionString):
    ctg,start,stop = re.search('(.*):(\d+)-(\d+)',regionString).groups()
    return ctg.strip(),int(start),int(stop)

def nameFunction(readname,start,stop):
    return f'{readname}/{start}_{stop}'

def getFlanks(ref,ctg,start,stop,flanksize=100,Lflank=None,Rflank=None):
    '''Extract flanking sequence from BED regions, given a pysam.FastaFile with .fai'''
    Lsize    = Lflank if Lflank else flanksize
    Rsize    = Rflank if Rflank else flanksize
    sequence = ref.fetch(ctg,start-Lsize,stop+Rsize)
    return [sequence[:Lsize],sequence[-Rsize:]]

def getFlankAligner(ref,ctg,start,stop,**kwargs):
    tmpRef = NamedTemporaryFile(mode='w',delete=False)
    for side,seq in zip(['L','R'],getFlanks(ref,ctg,start,stop,**kwargs)):
        tmpRef.write(f'>{"_".join([str(ctg),side])}\n{seq}\n')
    tmpRef.close()
    aligner = mp.Aligner(tmpRef.name,preset='sr')
    return aligner,tmpRef

def getSubSeq(seq,aln):
    pos = sorted([getattr(a,att) for a in aln for att in ['q_st','q_en']])[1:-1]
    return pos + [seq[slice(*pos)]]

def extractRepeat(sequence,aligner):
    aln = list(aligner.map(sequence))
    naln = len(aln)
    if naln == 2 and len(set(a.ctg for a in aln)) == 2:
        start,stop,seq = getSubSeq(sequence,aln)
        #seq = getSubSeq(sequence,aln)
        if aln[0].strand == -1:
            seq = rc(seq)
    else:
        seq = 'One Sided' if naln==1 else 'Poor/no Alignment'
        start,stop = None,None
    return start,stop,seq

_RC_MAP = dict(list(zip('-ACGNTacgt','-TGCNAtgca')))

def rc(seq):
    '''revcomp'''
    return "".join([_RC_MAP[c] for c in seq[::-1]])

def phred2prob(q):
    return 1-10**(-q/10)

def meanRQ(qual):
    return np.mean(list(map(phred2prob,map(int,qual))))

class Extract_Exception(Exception):
    pass
