#! /usr/bin/python

import sys,re
from collections import Counter
from pbcore.io import FastaReader,FastqReader

def main(parser):
    args = parser.parse_args()
    
    fx      = args.inFastx if args.inFastx else sys.stdin
    #sortedRecs    = sorted(pysam.FastxFile(fx),key=lambda r:-len(r.sequence))
    keyfunc = (lambda s:-len(s)) if args.reverse else len
    try:
        sortedRecs = sorted(FastqReader(fx),key=keyfunc)
    except ValueError:
        #this will fail if fasta is streamed 
        sortedRecs = sorted(FastaReader(fx),key=keyfunc)

    motifs    = args.motifs.split(',') 
    motifCols = args.sep.join(map('{{{}}}'.format,motifs))
    outFormat = '{{readName}}{sep}{cols}{sep}{{totalLength}}'.format(sep=args.sep,cols=motifCols)
    getCounts = countMotifs(motifs)

    oFile = open(args.out,'w') if args.out else sys.stdout
    #column names
    oFile.write(re.sub('{|}','',outFormat) + '\n')
    for rec in sortedRecs:
        counts = getCounts(rec.sequence)
        #fill for motifs not found
        counts.update({m:0 for m in motifs if not counts.has_key(m)})
        counts['readName']    = rec.name
        counts['totalLength'] = len(rec.sequence)
        oFile.write(outFormat.format(**counts) + '\n')
    oFile.close()
        
    return None

def countMotifs(motifs):
    '''returns a function that takes a sequence and returns a dict of counts'''
    patt = re.compile('(' + ')|('.join(motifs) + ')')
    def getCounts(seq):
        return Counter(m.group() for m in patt.finditer(seq))
    return getCounts

class CountMotifs_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='countMotifs.py', description='quick count of motifs per read from fastx of extracted repeat sequences')
    parser.add_argument('-i,--inFastx', dest='inFastx', type=str, default=None, 
                    help='Input Fastx file. Default stdin')
    parser.add_argument('-o,--out', dest='out', type=str, default=None,
                    help='Output csv. Default stdout')
    parser.add_argument('-m,--motifs', dest='motifs', type=str, default=None, required=True,
                    help='Search motifs, comma separated, most frequent first, e.g. \'CGG,AGG\'')
    parser.add_argument('-s,--sep', dest='sep', type=str, default=',',
                    help='Field separator.  Default \',\'')
    parser.add_argument('-r,--reverse', dest='reverse', action='store_true', default=False,
                    help='Sort largest first.  Default ascending order')

    try:
        main(parser)
    except CountMotifs_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1) 

