#! /usr/bin/python

from resources.utils import getFlankAligner,extractRepeat
import sys,re,pysam

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

