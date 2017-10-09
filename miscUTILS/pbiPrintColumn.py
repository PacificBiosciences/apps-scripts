#! /usr/bin/python

from pbcore.io import PacBioBamIndex
import sys

pbi=sys.argv[1]
col=sys.argv[2:]

print'\n'.join(map(lambda v:','.join(map(str,v)),zip(*[getattr(PacBioBamIndex(pbi),c) for c in col])))
