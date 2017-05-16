#!/usr/bin/env python

# Sarah B. Kingan
# 13 March 2017
# 
# Pacific Biosciences
# Applications Lab

###########################################################################################################

# import libraries
from argparse import ArgumentParser

###########################################################################################################

# define command line arguments, program description, and help
desc ='''Calculate mean QV for each contig.'''
parser = ArgumentParser(description=desc)
parser.add_argument("fastQ", help="FastQ file, no wrapping")
args = parser.parse_args()


###########################################################################################################

# get filenames
infile = args.fastQ


###########################################################################################################

# Iterative mean #
def imean(numbers):
    count = 0
    total = 0
    for num in numbers:
        count += 1
        total += num
    return float(total)/count

def wmean(means, lengths):
	m = 0
	t = 0
	for i in range(len(means)):
		m = m + (means[i] * lengths[i])
		t = t + lengths[i]
	return float(m)/t


# header
print "# @contigID"
print "#meanBaseQV\tContiglength"



l = 0
means = []
lengths = []
alist = [line.rstrip() for line in open(infile)]
for a in alist:
	l = l + 1
	m = l % 4
	if m == 1:
		print a
	elif m == 0:
		q = [ord(c)-33 for c in a]
		mq = imean(q)
		nq = len(a)
		print '\t'.join(str(o) for o in [mq, nq])
		means.append(mq)
		lengths.append(nq)
print '\t'.join(str(g) for g in ['global mean', wmean(means,lengths)])
