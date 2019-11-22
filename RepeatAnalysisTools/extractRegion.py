#! /usr/bin/python

import sys
from resources.extract import extractRegion,Extract_Exception

def main(parser):
    args = parser.parse_args()

    with (open(args.outFq,'w') if args.outFq else sys.stdout) as oFile:
        for name,seq,qual in extractRegion(args.inBAM,
                                           args.reference,
                                           region=args.region,
                                           flanksize=args.flanksize,
                                           revcomp=args.revcomp):
            oFile.write(fqRec(name,seq,qual))
    
    return None

def fqRec(name,seq,qual):
    return '@{name}\n{seq}\n+\n{qual}\n'.format(**locals())

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
    parser.add_argument('-r,--revcomp', dest='revcomp', action='store_true', default=False,
                    help='Rev-comp extracted region.  Default Reference Direction')

    try:
        main(parser)
    except Extract_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1) 

