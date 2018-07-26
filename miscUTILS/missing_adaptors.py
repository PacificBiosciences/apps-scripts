#!/usr/local/env python
"""Calculate alignment of sequences flanking adaptor and plot distribution.
"""

import os.path
import sys
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
import numpy as np
import pandas as pd

from Bio import pairwise2, SeqIO


FLANK_LEN = 100
MATCH = 5
MISMATCH = -4
GAP_OPEN = -10
GAP_EXTEND = -0.5


def calc_identity(record):
    """Calculate the percent identity of adaptor flanking sequences.

    refer to http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    """
    leftflank = record.seq[0:FLANK_LEN]
    rightflank = record.seq[-FLANK_LEN:].reverse_complement()
    alignment = pairwise2.align.globalms(leftflank,
                                         rightflank,
                                         MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND,
                                         score_only=True)
    return alignment, (float(alignment) / (MATCH * FLANK_LEN)) * 100


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fasta', help='input fasta', type=str)
    args = parser.parse_args(arguments)

    outfig = os.path.splitext(args.fasta)[0] + '.pdf'
    outcsv = os.path.splitext(args.fasta)[0] + '.csv'

    data = defaultdict(list)
    for record in SeqIO.parse(args.fasta, 'fasta'):
        data['name'].append(record.id)
        data['zmw'].append(record.id.split('/')[1])
        data['length'].append(len(record))
        score, percent = calc_identity(record)
        data['score'].append(score)
        data['percent'].append(percent)
    df = pd.DataFrame(data=data)
    df = df[['name', 'zmw', 'length', 'score', 'percent']]

    # plot distributions
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    df.percent.hist(ax=ax, bins=np.linspace(df.percent.min(), df.percent.max(), 50), log=True)
    ax.set_xlabel('percent alignment')
    ax.set_ylabel('log count')
    plt.tight_layout()
    plt.savefig(outfig)

    # save table
    df.to_csv(outcsv, sep='\t', index=False)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
