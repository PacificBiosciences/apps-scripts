#! /usr/bin/python

import sys,re,pysam
import mappy as mp
from tempfile import NamedTemporaryFile

ALIGNFILTER=0x900

def main(parser):
    args = parser.parse_args()

    ref = pysam.FastaFile(args.reference)
    bam = pysam.AlignmentFile(args.inBAM)
    try:
        ctg,start,stop = getCoordinates(args.region)
    except AttributeError:
        #catch when the region format doesn't match
        raise ExtractRegion_Exception('Invalid region format %s. Correct \'[chr]:[start]-[stop]\'' % args.region)

    aligner,tmp = getFlankAligner(ref,ctg,start-1,stop,flanksize=args.flanksize) 
    
    with (open(args.outFq,'w') if args.outFq else sys.stdout) as oFile:
        for rec in bam.fetch(ctg,start,stop):
            if (rec.flag & ALIGNFILTER):
                continue
            rStart,rStop,subseq = extractRepeat(rec.query_sequence,aligner)
            if rStart:
                name = nameFunction(rec.query_name,rStart,rStop)
                qual = ''.join([chr(q+33) for q in rec.query_qualities[rStart:rStop]])
                oFile.write(fqRec(name,subseq,qual))

def getCoordinates(regionString):
    ctg,start,stop = re.search('(.*):(\d+)-(\d+)',regionString).groups()
    return ctg.strip(),int(start),int(stop)

def nameFunction(readname,start,stop):
    return '{readname}/{start}_{stop}'.format(**locals())

def fqRec(name,seq,qual):
    return '@{name}\n{seq}\n+\n{qual}\n'.format(**locals())

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

class ExtractRegion_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='extractRegion.py', description='extract target region from aligned BAMS using region flank alignments. Output format is fastq')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='input BAM of CCS alignments')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta used for mapping BAM.  Must have .fai index.')
    parser.add_argument('region', metavar='region', type=str, 
                    help='Target region, format \'[chr]:[start]-[stop]\'.  Example \'4:3076604-3076660\'')
    parser.add_argument('-o,--outFq', dest='outFq', type=str, default=None,
                    help='Output fastq file.  Default stdout')
    parser.add_argument('-f,--flanksize', dest='flanksize', type=int, default=100,
                    help='Size of flanking sequence mapped for extracting repeat region.  Default 100')

    try:
        main(parser)
    except ExtractRegion_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1) 

