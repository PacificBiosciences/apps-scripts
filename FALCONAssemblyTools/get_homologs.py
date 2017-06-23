#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple python script to identify homologous contigs in FALCON primary contig output
Greg Concepcion gconcepcion@pacificbiosciences.com
"""

import os
import csv
import sys
import glob
import argparse
import subprocess
import logging
from collections import Counter

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def explode_fasta(fasta):
    """split input fasta into individuals"""
    in_file = False
    outdir = "fastas"

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    with open(fasta, 'r') as infile:

        for line in infile:
            if line.startswith(">"):
                if in_file:
                    outfile.close()

                fname = line.rstrip().partition(">")[2].rstrip('|arrow')
                fname = "{s}.fasta".format(s=fname)
                fout = os.path.join(outdir, fname)
                outfile = open(fout, "w")
                in_file = True

                outfile.write(line)
            elif in_file:
                outfile.write(line)
            else:
                log.debug("Line %r, but no previous > found ")

    return glob.glob('fastas/*')


def run_nucmer(reference, queries, threads):
    """run nucmer for each reference against all queries"""
    refname = os.path.basename(reference.rstrip(".fasta"))

    log.info("Searching all queries for alignments with reference %s", refname)

    if not os.path.exists('deltas'):
        os.mkdir('deltas')
    prefix = "deltas/{f}".format(f=refname)

    cmd = ["nucmer", "--maxmatch", "-l", "100", "-c", "500",
           "-t", str(threads), "-p", prefix, reference, queries]
    log.debug(cmd)
    run(cmd)

    deltafile = "{p}.delta".format(p=prefix)
    return deltafile


def run_show_coords(deltafile):
    """run show-coords for each reference"""
    log.debug("Converting delta to coords")

    cmd = ["show-coords", "-HT", deltafile]
    log.debug(cmd)
    stdout = run(cmd)

    coords_file = [line for line in stdout.split(os.linesep)]

    return coords_file


def run_delta_filter(deltafile):
    """Generate filtered delta file for mummerplot"""

    log.debug("filtering delta for plotting")
    new_delta = "{d}_filtered.delta".format(d=deltafile.rstrip('.delta'))
    cmd = ["delta-filter", "-g", deltafile]

    stdout = run(cmd)

    with open(new_delta, 'w') as output:
        for line in stdout:
            output.write(line)

    return new_delta


def run(cmd):
    """Run shell process"""

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stderr:
        log.debug(stderr)

    return stdout.rstrip()


def process_coords(coordsfile, delta, length_dict):
    """Process *.coords output for significant matches"""

    log.debug("Processing coords file")
    reference = os.path.basename(delta.rstrip('.delta'))

    coords = [tuple(i.split()) for i in coordsfile]
    count = Counter(i[8] for i in coords)

    filtered = [i for i, c in count.iteritems() if c > 3]

    query_dict = {}
    self = None
    for query in filtered:
        query_hits = []

        for line in coords:
            if line[7] == line[8]:
                self = line[7]
                #pass
            elif line[8] == query:
                query_hits.append(line)
        query_dict[query] = query_hits
        query_dict.pop(self, None)

    new_list = []

    for key, value in query_dict.iteritems():
        startends = []
        total_bp = sum([int(i[4]) for i in value])
        percent_ref = round(total_bp / float(length_dict[reference]), 4)

        for hit in value:
            startends.append((hit[0], hit[1]))

        if len(merge(startends)) / float(len(startends)) > 0.75:
            if percent_ref >= 0.03:
                ratio = len(merge(startends)) / float(len(startends))
                new_list.append((key, total_bp, percent_ref, ratio))

    log.info("%s shares homology with %s", reference, ",".join([i[0] for i in new_list]))
    qfile = write_qfile(reference, new_list, length_dict)

    return qfile


def write_qfile(reference, contigs, length_dict):
    """Write qfile of homologous IDs for mummerplot"""

    qfiledir = os.path.join(os.getcwd(), 'qfiles')
    if not os.path.exists(qfiledir):
        os.mkdir(qfiledir)
    qfile = os.path.join(qfiledir, "{r}.qfile".format(
        r=os.path.basename(reference.rstrip('.fasta'))))

    with open(qfile, 'w') as fout:

        csv_out = csv.writer(fout, delimiter=' ', lineterminator="\n")

        for contig in contigs:
            name = contig[0].rstrip('|arrow')
            length = length_dict[name]
            csv_out.writerow((contig[0], length, '+'))

    return qfile


def get_length_dict(fastas):
    """Generate length dictionary for all contigs"""

    length_dict = {}

    for fasta in fastas:
        name = os.path.basename(fasta.rstrip('.fasta'))
        length = get_length(fasta)
        length_dict[name] = length

    return length_dict


def get_length(fastafile):
    """Get length of a fasta sequence"""

    length = int()
    with open(fastafile, 'r') as fin:
        for line in fin:

            if line.startswith(">"):
                pass
            else:
                length += len(line.strip())
        return length


def merge(qhits):
    """Merge overlapping hits together to get a contiguous interval"""

    intervals = [(int(i[0]), int(i[1])) for i in qhits]

    if not intervals:
        return []
    data = []
    for interval in intervals:
        data.append((interval[0], 0))
        data.append((interval[1], 1))
    data.sort()
    merged = []
    stack = [data[0]]
    for i in xrange(1, len(data)):
        datum = data[i]
        if datum[1] == 0:
            # this is a lower bound, push this onto the stack
            stack.append(datum)
        elif datum[1] == 1:
            if stack:
                start = stack.pop()
            if len(stack) == 0:
                # we have found our merged interval
                merged.append((start[0], datum[0]))
    return merged


def get_mummerplot_cmd(delta, qfile):
    """Generate mummerplot command for each reference"""
    prefix = os.path.basename(delta).rstrip('.delta')
    new_delta = run_delta_filter(delta)

    cmd = "mummerplot -layout -Q {q} -postscript -p {p} {d}".format(q=qfile,
                                                                    d=new_delta,
                                                                    p=prefix)
    return cmd


def write_bed(reference, query, qhits):
    """Write *.bed annotation file"""

    ref = reference.split('|')[0]
    que = query.split('|')[0]
    outdir = os.path.join(os.getcwd(), 'bedfiles')

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    bedname = "{r}_{q}.bed".format(r=ref, q=que)
    bedfile = os.path.join(outdir, bedname)

    with open(bedfile, 'w') as out:
        csv_out = csv.writer(out, delimiter='\t')

        for row in qhits:
            rstart, rend, _, _, _, _, _ = row
            csv_out.writerow((reference, int(rstart), int(rend), query))

    return bedfile


def get_parser():
    """Return an argparse instance"""

    __version__ = 0.1
    parser = argparse.ArgumentParser(version=__version__)
    parser.add_argument("infile", type=str)
    parser.add_argument("--nproc", type=int, default=8)
    parser.add_argument('--debug', action='store_true',
                        help="Print debug logging to stdout")

    return parser.parse_args()


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s ' \
                            '[%(name)s %(funcName)s %(lineno)d] ' \
                            '%(message)s'):
    """Core Util to setup log handler"""
    alog.setLevel(logging.DEBUG)
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)
    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    if log_filter:
        handler.addFilter(log_filter)
    alog.addHandler(handler)


def main():
    """Main run loop"""

    args = get_parser()
    infile = args.infile
    threads = args.nproc
    debug = args.debug
    if debug:
        setup_log(log, file_name='log.out', level=logging.DEBUG)
    else:
        setup_log(log, file_name='log.out', level=logging.INFO)

    if infile.endswith(('.fasta', '.fa')):
        fastas = explode_fasta(infile)
    else:
        log.info("Please provide FASTA as your input file")

    length_dict = get_length_dict(fastas)
    total_seqs = len(length_dict.keys())
    length_sum = sum(length_dict.values())
    log.info("Beginning homology search in %s", infile)
    log.info("Total Contigs: %d", total_seqs)
    log.info("Total Bp: %d", length_sum)

    plot_out = 'plots.sh'

    with open(plot_out, 'w') as plot_out:
        for fasta in sorted(fastas):

            delta = run_nucmer(fasta, infile, threads)
            coordsfile = run_show_coords(delta)

            qfile = process_coords(coordsfile, delta, length_dict)
            plot_cmd = get_mummerplot_cmd(delta, qfile)

            plot_out.write(plot_cmd + '\n')


if __name__ == "__main__":
    sys.exit(main())
