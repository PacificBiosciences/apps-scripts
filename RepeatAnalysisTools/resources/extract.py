import os,re,pysam
import mappy as mp
from tempfile import NamedTemporaryFile

ALIGNFILTER=0x900

def extractRegion(inBAM,reference,region=None,ctg=None,start=None,stop=None,flanksize=100,revcomp=False):
    ref = pysam.FastaFile(reference)
    bam = pysam.AlignmentFile(inBAM)
    if region:
        try:
            ctg,start,stop = getCoordinates(region)
        except AttributeError:
            #catch when the region format doesn't match
            raise Extract_Exception('Invalid region format %s. Correct \'[chr]:[start]-[stop]\'' % region)
    elif ctg==None or start==None or stop==None:
            #catch when missing coord and no region passed
            raise Extract_Exception('Must pass either valid region string or all of ctg,start,stop')

    aligner,tmp = getFlankAligner(ref,ctg,start-1,stop,flanksize=flanksize) 
    
    try:
        for rec in bam.fetch(ctg,start,stop):
            if (rec.flag & ALIGNFILTER):
                continue
            rStart,rStop,subseq = extractRepeat(rec.query_sequence,aligner)
            if rStart:
                name = nameFunction(rec.query_name,rStart,rStop)
                qual = ''.join([chr(q+33) for q in rec.query_qualities[rStart:rStop]])
                if revcomp:
                    subseq = rc(subseq)
                    qual   = qual[::-1]
                yield name,subseq,qual
    finally:
        os.remove(tmp.name)

def getCoordinates(regionString):
    ctg,start,stop = re.search('(.*):(\d+)-(\d+)',regionString).groups()
    return ctg.strip(),int(start),int(stop)

def nameFunction(readname,start,stop):
    return '{readname}/{start}_{stop}'.format(**locals())

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
    aligner = mp.Aligner(tmpRef.name,preset='sr')
    return aligner,tmpRef

def getSubSeq(seq,aln):
    pos = sorted([getattr(a,att) for a in aln for att in ['q_st','q_en']])[1:-1]
    return pos + [seq[slice(*pos)]]

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

_RC_MAP = dict(zip('-ACGNTacgt','-TGCNAtgca'))

def rc(seq):
    '''revcomp'''
    return "".join([_RC_MAP[c] for c in seq[::-1]])

class Extract_Exception(Exception):
    pass
