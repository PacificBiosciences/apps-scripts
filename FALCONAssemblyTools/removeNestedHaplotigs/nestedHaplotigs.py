#!/usr/bin/env python

# Sarah B. Kingan
# 7 April 2017
# 
# Pacific Biosciences
# Applications Lab

###########################################################################################################

# import libraries
from argparse import ArgumentParser
import csv
import numpy as np
import subprocess
import string

###########################################################################################################

# define command line arguments, program description, and help
desc ='''Identify nested haplotigs.'''
parser = ArgumentParser(description=desc)
parser.add_argument("fileList", help="list of coords files")
args = parser.parse_args()


###########################################################################################################

# get filenames
inFiles = args.fileList

###########################################################################################################

# define functions

def file2list(textFile):
	return [line.rstrip('\n') for line in open(textFile)]

def file_len(fname):
	#http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
	p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	result, err = p.communicate()
	if p.returncode != 0:
		raise IOError(err)
	return int(result.strip().split()[0])

def bed(file):
	if file_len(file) >= 5:
		output=['primaryContig','pstart', 'pend','haplotig']
		d=np.loadtxt(open(file, "rb"), skiprows=4, dtype="str", ndmin=2)
		output[3] = d[0,10] # haplotig
		output[0] = d[0,9] #primary contig
		p = np.array(d[:,[0,1]],dtype=int)
		output[1] = np.amin(p) # pstart
		output[2] = np.amax(p) # pstop
		return output
	
def nested(a,b):
	# a nested inside b
	nested = 0
	if a[1] < b[2] and a[1] > b[1] and a[2] > b[1] and a[2] < b[2]:
		nested = 1
	return nested


coordsFileList = file2list(inFiles)
nestedHaplotigs = [0] * len(coordsFileList)
for i in range(len(coordsFileList)):
	for j in range(i):
		bedi = bed(coordsFileList[i])
		bedj = bed(coordsFileList[j])
		if nested(bedi,bedj):
			nestedHaplotigs[i] = 1
		if nested(bedj,bedi):
			nestedHaplotigs[j] = 1

# ready to write to file
nested = open('nested_haplotigs.txt', 'a')
retained = open('retained_contigs_haplotigs.txt', 'a')

for k in range(len(coordsFileList)):
	haplotig = str(coordsFileList[k])
	haplotig = haplotig.replace(".coords","")
	if nestedHaplotigs[k] == 1:
		nested.write(str(haplotig) + "\n")
	else:
		retained.write(str(haplotig) + "\n")

	
nested.close()
retained.close()
