import pysam,sys
import numpy as np
from collections import Counter
from resources.extract import extractRegion,fqRec,rc

MINCLUSTERSIZE=5

def main(parser):
    args = parser.parse_args()

    if args.inBAM and args.inFastq:
        raise SampleReads_Exception('Only one input, either -b or -q') 

    if args.inBAM:
        bam = pysam.AlignmentFile(args.inBAM,check_sq=False)
        names = getReadNamesBam(bam,HP=args.haplotag)
        if args.region:
            if not args.reference:
                raise SampleReads_Exception('Must pass reference for region extraction')

            recGen = extractRegion(args.inBAM,
                                   args.reference,
                                   region=args.region,
                                   flanksize=args.flanksize,
                                   revcomp=args.revcomp)
        else:
            recGen = recIterBam(bam,revcomp=args.revcomp)
    elif args.inFastq:    
        names  = getReadNamesFq(args.inFastq)
        recGen = recIterFq(args.inFastq,revcomp=args.revcomp) 

    else:
        raise SampleReads_Exception('Must have input, either -b or -q') 
        

    if not len(names):
        raise SampleReads_Exception('No reads returned')
    if len(names)<MINCLUSTERSIZE:
        raise SampleReads_Exception('Fewer than %i reads returned' % MINCLUSTERSIZE)


    np.random.seed(args.seed)
    size     = args.nReads if args.nReads else len(names)
    selected = Counter(np.random.choice(names,size=size,replace=args.replace))


    nrecs = 0
    with (open(args.out,'w') if args.out else sys.stdout) as oFile:
        for name,seq,qual in recGen:
            cname = clipReadName(name)
            if cname in selected:
                rec = fqRec(cname,seq,qual)
                nrecs += 1
                for _ in range(selected[cname]):
                    oFile.write(rec)
            if nrecs == size:
                break
    return None

def recIterBam(bam,revcomp=False):
    for rec in bam:
        if rec.flag & 0x900:
            continue
        seq  = rec.query_sequence
        qual = ''.join([chr(q+33) for q in rec.query_qualities])
        if revcomp:
            seq  = rc(seq)
            qual = qual[::-1]
        yield rec.query_name,seq,qual

def recIterFq(fastq,revcomp=False):
    trs = rc if revcomp else (lambda x:x)
    trq = (lambda v:v[::-1]) if revcomp else (lambda x:x)
    for rec in pysam.FastxFile(fastq):
        yield rec.name,trs(rec.sequence),trq(rec.quality)

def getReadNamesBam(bam,HP=None):
    crit = (lambda rec: True) if HP is None else (lambda rec: rec.get_tag('HP')==HP)
    try:
        return sorted(set(rec.query_name for rec in bam if crit(rec)))
    except KeyError:
        raise SampleReads_Exception('No HP tag in BAM')
    finally:
        bam.reset()

def getReadNamesFq(fastq):
    return [clipReadName(rec.name) for rec in pysam.FastxFile(fastq)]

def clipReadName(name,nFields=3):
    return '/'.join(name.split('/')[:nFields])

class SampleReads_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='sampleBam.py', description='export a random sample of reads in fastq format')
    parser.add_argument('-b','--inBAM', dest='inBAM', type=str, default=None,
                    help='BAM containing reads to sample.  Default None')
    parser.add_argument('-q','--inFastq', dest='inFastq', type=str, default=None,
                    help='fastq containing reads to sample.  No region or HP filtering, use all reads.  Default None')
    parser.add_argument('-n','--nReads', dest='nReads', type=int, default=0,
                    help='Reads to sample.  To resample for bootstrapping, use default. Default 0 (all reads)')
    parser.add_argument('--reg', dest='region', type=str, default=None,
                    help='Target region to extract, format \'[chr]:[start]-[stop]\'.  Example \'4:3076604-3076660\'. Default None.')
    parser.add_argument('--ref', dest='reference', type=str, default=None,
                    help='Reference fasta used for mapping BAM if extracting region.  Must have .fai index. Default None')
    parser.add_argument('-f','--flanksize', dest='flanksize', type=int, default=100,
                    help='Size of flanking sequence mapped for extracting repeat region.  Default 100')
    parser.add_argument('--rc', dest='revcomp', action='store_true', default=False,
                    help='Rev-comp extracted region.  Default Reference Direction')    
    parser.add_argument('-H','--haplotag', dest='haplotag', type=int, default=None,
                    help='Sample from one HP tag value. Default None (all reads)')
    parser.add_argument('-s','--seed', dest='seed', type=int, default=17,
                    help='Random seed. Default 17')
    parser.add_argument('-o','--out', dest='out', type=str, default=None,
                    help='Output file. Default stdout')
    parser.add_argument('-r','--replace', dest='replace', action='store_true', default=False,
                    help='Sample with replacement. default False')

    try:
        main(parser)
    except SampleReads_Exception as e:
        print('ERROR: %s' % e)
        sys.exit(1)
