#! /usr/bin/python

import sys,re
from collections import Counter,OrderedDict
from pbcore.io import FastaReader,FastqReader
from resources.utils import countMotifs

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
    if len(sortedRecs) == 0:
        #weird single fx rec with three lines fails
        sortedRecs = sorted(FastaReader(fx),key=keyfunc)

    motifs    = args.motifs.split(',')
    motifCols = args.sep.join(map('{{{}}}'.format,motifs))
    if args.blockCounts:
        outFormat = f'{{readName}}{args.sep}{motifCols}{args.sep}{{blockCount}}{args.sep}{{totalLength}}'
    else:
        outFormat = f'{{readName}}{args.sep}{motifCols}{args.sep}{{totalLength}}'
    getCounts = countMotifs(motifs,
                            lengthField='totalLength',
                            blocks='blockCount' if args.blockCounts else False,
                            collapseHP=args.collapseHP)

    oFile = open(args.out,'w') if args.out else sys.stdout
    #column names
    oFile.write(re.sub('{|}','',outFormat) + '\n')
    for rec in sortedRecs:
        counts = getCounts(rec.sequence)
        counts['readName']    = rec.name
        oFile.write(outFormat.format(**counts) + '\n')
    oFile.close()
        
    return None

_PATT = re.compile(r'([ATGC])\1+')
def hpCollapse(seq):
    return _PATT.sub(r'\1',seq)

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
    parser.add_argument('-b,--blockCounts', dest='blockCounts', action='store_true', default=False,
                    help='Count primary motif blocks with secondary as interruption.  Default false')
    parser.add_argument('-r,--reverse', dest='reverse', action='store_true', default=False,
                    help='Sort largest first.  Default ascending order')
    parser.add_argument('-c,--collapseHP', dest='collapseHP', action='store_true', default=False,
                    help='Count motifs after HP collapsing both motif and sequence (experimental)')

    try:
        main(parser)
    except CountMotifs_Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1) 

