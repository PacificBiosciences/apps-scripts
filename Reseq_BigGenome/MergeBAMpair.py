#!/usr/bin/env python

# Sarah B. Kingan
# 20 October 2017
# updated 21 November 2017
# 
# Pacific Biosciences
# BFX Applications

###########################################################################################################

# import libraries
import os
import pysam
import csv
from argparse import ArgumentParser
from collections import defaultdict

###########################################################################################################

def make_readID_hash(file):
        d = defaultdict(dict)
        with open(file, mode='r') as infile:
                reader = csv.reader(infile, delimiter='/')
                for rows in reader:
	                        d[rows[0]][rows[1]] = {'AS1': 'AS1', 'AS2': 'AS2', 'N1': 'N1', 'N2': 'N2'}
        return d

def get_movie(readID):
	myList = readID.split("/")
	return str(myList[0])

def get_zmw(readID):
	myList = readID.split("/")
	return str(myList[1])

###########################################################################################################

# define command line arguments, program description, and help
desc ='''Input BAM files and list of "shared reads". Output BAM files with 
deduplicated reads.'''
parser = ArgumentParser(description=desc)
parser.add_argument("BAM1", help="BAM file for contigset 1")
parser.add_argument("BAM2", help="BAM file for contigset 2")
args = parser.parse_args()

# get infiles and dir
cwd = os.getcwd()
sharedReads = str(cwd + '/sharedReads.txt')
dest = str(cwd + '/newBAMs/contigSet')
BAM1 = args.BAM1
BAM2 = args.BAM2

# check infiles are same chunk (same subread set)
if BAM1.split("/")[-2].split("-")[-1] != BAM2.split("/")[-2].split("-")[-1]:
	print "ERROR: BAMs have different subreadsets"
	sys.exit()

# subreadset chunk ID
fileID = BAM1.split("/")[-2].split("-")[1]

# out BAM filenames
out1 = ".".join(str(o) for o in [dest,"1",fileID,"bam"])
out2 = ".".join(str(o) for o in [dest,"2",fileID,"bam"])

# get shared read IDs into dictionary
sharedReadDict = make_readID_hash(sharedReads)

# initialize pysam objects for inputs and output
inbam1 = pysam.AlignmentFile(BAM1,'rb')
outbam1 = pysam.AlignmentFile(out1, "wb", template=inbam1)
inbam2 = pysam.AlignmentFile(BAM2,'rb')
outbam2 = pysam.AlignmentFile(out2, "wb", template=inbam2)


# loop through contiset1 BAM
for read in inbam1.fetch():
# populate shared read dict with alignment score data

	movie = get_movie(read.query_name)
	zmw = get_zmw(read.query_name)
	AS = int(read.get_tag('AS'))

	if movie in sharedReadDict and zmw in sharedReadDict[movie]:

		
		# initialize new aln score and count
		if sharedReadDict[movie][zmw]['AS1'] == 'AS1':
			sharedReadDict[movie][zmw]['AS1'] = AS
			sharedReadDict[movie][zmw]['N1'] = int(1)
		# update best aln score and count
		else:
			if AS < sharedReadDict[movie][zmw]['AS1']:
				sharedReadDict[movie][zmw]['AS1'] = AS
			sharedReadDict[movie][zmw]['N1'] += 1

# or print non-shared reads to file
	else:
        	outbam1.write(read)


# loop through contiset2 BAM
for read in inbam2.fetch():
# populate shared read dict with alignment score data

        movie = get_movie(read.query_name)
        zmw = get_zmw(read.query_name)
        AS = int(read.get_tag('AS'))

        if movie in sharedReadDict and zmw in sharedReadDict[movie]:


                # initialize new aln score and count
                if sharedReadDict[movie][zmw]['AS2'] == 'AS2':
                        sharedReadDict[movie][zmw]['AS2'] = AS
                        sharedReadDict[movie][zmw]['N2'] = int(1)
                # update best aln score and count
                else:   
                        if AS < sharedReadDict[movie][zmw]['AS2']:
                                sharedReadDict[movie][zmw]['AS2'] = AS
                        sharedReadDict[movie][zmw]['N2'] += 1

# or print non-shared reads to file
        else:
                outbam2.write(read)




# initialize pysam objects for inputs and output for pass 3
inbam3 = pysam.AlignmentFile(BAM1,'rb')
out3 = ".".join(str(o) for o in [dest,"1s",fileID,"bam"])
outbam3 = pysam.AlignmentFile(out3, "wb", template=inbam3)

for read in inbam1.fetch():
# check if contigset has best aln score

        movie = get_movie(read.query_name)
        zmw = get_zmw(read.query_name)

        if movie in sharedReadDict and zmw in sharedReadDict[movie]:
                if sharedReadDict[movie][zmw]['AS1'] <  sharedReadDict[movie][zmw]['AS2']:
			outbam3.write(read)


# initialize pysam objects for inputs and output for pass 4
inbam4 = pysam.AlignmentFile(BAM2,'rb')
out4 = ".".join(str(o) for o in [dest,"2s",fileID,"bam"])
outbam4 = pysam.AlignmentFile(out4, "wb", template=inbam4)

for read in inbam2.fetch():
# check if contigset has best aln score

        movie = get_movie(read.query_name)
        zmw = get_zmw(read.query_name)

        if movie in sharedReadDict and zmw in sharedReadDict[movie]:
                if sharedReadDict[movie][zmw]['AS2'] <  sharedReadDict[movie][zmw]['AS1']:
                        outbam4.write(read)

# print info about shared reads
outdict = ".".join(["SharedReadDict",fileID,"txt"])
f = open(outdict, "w")
for movie in sharedReadDict:
	for zmw in sharedReadDict[movie]:
		if sharedReadDict[movie][zmw]['N1'] != 'N1' or sharedReadDict[movie][zmw]['N2'] != 'N2':
			output = [movie, zmw, sharedReadDict[movie][zmw]['AS1'], sharedReadDict[movie][zmw]['N1'], sharedReadDict[movie][zmw]['AS2'], sharedReadDict[movie][zmw]['N2']]
			f.write("\t".join(str(o) for o in output))
			f.write("\n")
