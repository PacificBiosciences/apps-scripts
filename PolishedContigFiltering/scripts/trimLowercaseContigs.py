#!/usr/bin/env python

# This tool requires the pacbio pbcore library
# Filter out contigs that are made up of a defined percentage of lowercase charagters, for Pacbio HGAP.4 
# assemblies this is an indication that an arrow consensus was not calculated for that particular base.
# Written by Richard Hall

import sys
import os
import csv
import re
import subprocess
from shutil import copy2
from collections import namedtuple
import argparse

from pbcore.io.FastqIO import FastqReader, FastqWriter, FastqRecord
from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def main(argv):
        desc = 'A tool to trim quiver results for contigs majority lowercase'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('inputFile',  help='input sequence')
        parser.add_argument('outputFile', help='output fasta')
        parser.add_argument('--filt', default=0.5,
                            dest='filt', type=float,
                            help='proportion of lowercase bases a contig can have before being filtered out')
        args = parser.parse_args()

        writer = FastaWriter(args.outputFile)

        for record in FastaReader( args.inputFile ):
            upper_output = []
            upper_indx = []
            lower = float(sum(1 for c in record.sequence if c.islower()))
            pro = lower / float(len(record.sequence))
            print pro
            if pro < args.filt:
                writer.writeRecord(record)

if __name__ == "__main__":
   main(sys.argv[1:])
