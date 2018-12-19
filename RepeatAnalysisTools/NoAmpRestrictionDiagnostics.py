#! /home/UNIXHOME/jharting/anaconda2/bin/python

import numpy as np
import pandas as pd
import re,os,sys
import pysam
from pbcore.io import FastaReader
from resources.utils import rc,\
                            readBED,\
                            getHn

ALIGNFILTER=0x900
FLOATFORMAT='%.2f'
PCTFUNC='{:,.2%}'.format
NOADAPTER='NoAd'

def main(parser):
    args = parser.parse_args()

    if args.subreadsBAM or args.adapterFasta:
        if not args.subreadsBAM and args.adapterFasta:
            raise NoAmpRestrictionDiagnostics_Exception, \
                  'Must use both -s and -a if classifying adapters'
        #parse adapter(ad) tags
        print 'Parsing Adapter information'
        bam      = pysam.AlignmentFile(args.subreadsBAM,
                                       check_sq=False)
        adapters = pd.DataFrame([{'holeNumber':getHn(rec.query_name),
                                  'ad'        :rec.get_tag('ad')}
                         for rec in bam 
                         if rec.has_tag('ad')])
        adParser = adapterParser(args.adapterFasta,
                                 minAlnScore=args.minAdapterScore,
                                 minFlankScore=args.minFlankScore)
        splitAds = adapters.ad.str.split(';',expand=True)

        adapters['ad_parsed'] = splitAds.apply(adParser)\
                                        .apply(lambda v: tuple(sorted(v)),axis=1)

        #classify parsed adapters
        adClassifier = classifyZmw(args.adapterFasta)
        holeClasses  = adapters.groupby('holeNumber')\
                               .ad_parsed.apply(adClassifier)\
                               .str.join(';')
        #build summary table (write later)
        noAdaps  = ';'.join([NOADAPTER,NOADAPTER])
        cnts     = holeClasses.value_counts()\
                              .append(pd.Series({noAdaps:bam.unmapped - len(adapters)}))\
                              .rename('ZMWcount')\
                              .sort_values(ascending=False)
        frac     = (1.*cnts/cnts.sum()).rename('TotalFraction')\
                                       .map(PCTFUNC)
        adreport = pd.DataFrame([cnts,frac]).T
        adreport.index.name = 'adapters'
    else:
        print 'WARNING: No adapter information will be generated'

    print 'Tabulating Restriction sites in reads'
    aBam    = pysam.AlignmentFile(args.inBAM)
    targets = readBED(args.inBED,
                      names=['ctg','start','end','name'])
    targets.set_index('name',inplace=True)
    ref     = pysam.FastaFile(args.reference)
    reTable = readRestrictionTsv(args.reTable)
    #for each read intersecting a target, get holenumber 
    #and ref start/stop of primary alignment
    onTarget = pd.DataFrame([{'holeNumber': getHn(rec.query_name),
                              'target'    : name,
                              'ctg'       : target.ctg,
                              'rStart'    : rec.reference_start,
                              'rEnd'      : rec.reference_end}
                             for name,target in targets.iterrows()
                             for rec in aBam.fetch(target.ctg,target.start,target.end)
                             if not (rec.flag & ALIGNFILTER)])
    #the part of the reference sequence in the region of each target that is
    #covered by at least one read (plus some extra offset/buffer)
    #we will search the covered regions of the reference for ny cutsite listed in the reTable
    offset           = reTable.cutsite.str.len().max()
    spannedReference = onTarget.groupby(['target','ctg'],as_index=False)\
                               .agg({'rStart': lambda s:max(0,s.min()-offset),
                                     'rEnd'  : lambda e:e.max()+offset})
    spannedReference['sequence'] = spannedReference.apply(lambda t: ref.fetch(t.ctg,
                                                                              t.rStart,
                                                                              t.rEnd),
                                                          axis=1)
    #initialize empty values
    #need to do this outside/first because 
    #some rows might be duplicates when
    #cutsites have degenerate bases
    #e.g. ACC1
    for enz in reTable.enzName:
        spannedReference[enz] = [[]]*len(spannedReference)
    
    #now we update each column with cutsite positions
    for i,rEnz in reTable.iterrows():
        #fwd
        spannedReference[rEnz.enzName] += spannedReference.sequence.map(getREpos(rEnz.cutsite))
        #revcomp
        spannedReference[rEnz.enzName] += spannedReference.sequence.map(getREpos(rc(rEnz.cutsite)))
    
    #after accumulating, return only unique locations
    #and drop enzyme cols with no found cutsites
    foundCutsites = set()
    for enz in reTable.enzName:
        if enz not in spannedReference.columns:
            #skip degenerate cutsite enzymes already dropped
            continue
        elif (spannedReference[enz].str.len() == 0).all():
            #drop if no cutsites found
            spannedReference.drop(enz,axis=1,inplace=True)
        else:
            #unique sites only
            spannedReference[enz] = spannedReference[enz].map(set)
            #add enz to found list
            foundCutsites.add(enz)
        
    #walk through found sites and count up occurrances in on-target reads
    #left end (lend) ; internal (int) ; right end (rend)
    suff = {'lend':'_lend','int':'_int','rend':'_rend'}
    for enz in foundCutsites:
        for s in suff.values():
            #initialize empty column counts
            onTarget[enz+s] = 0
        for target,row in spannedReference.iterrows():
            for start,stop in row[enz]:
                isOnTargetCtg = onTarget.ctg == row.ctg 
                lend = isOnTargetCtg & (onTarget.rStart >= start + row.rStart)\
                                     & (onTarget.rStart <= stop  + row.rStart)
                rend = isOnTargetCtg & (onTarget.rEnd   >= start + row.rStart)\
                                     & (onTarget.rEnd   <= stop  + row.rStart)
                mid  = isOnTargetCtg & (onTarget.rStart <  start + row.rStart)\
                                     & (onTarget.rEnd   >  stop  + row.rStart)
                #ad hits to table
                onTarget.loc[lend,enz+suff['lend']] += 1
                onTarget.loc[rend,enz+suff['rend']] += 1
                onTarget.loc[mid,enz+suff['int']]   += 1
    #drop columns (e.g. "EcoR1_lend") with no hits
    onTarget = onTarget.loc[:,~(onTarget==0).all()]
    #remaining enzyme columns are all to right of 'target'
    enzCols  = list(onTarget.columns[onTarget.columns.get_loc('target') + 1:])
    grpCols = ['target']
    if args.adapterFasta:
        #include adapter information from subread data
        onTarget['adapters'] = onTarget.holeNumber.map(holeClasses.to_dict())
        grpCols += ['adapters']
        #also write adapter report here, after adding CCS info
        aBam.reset()
        adreport['mappedCCS'] = holeClasses[[getHn(rec.query_name) for rec in aBam
                                             if not (rec.flag & ALIGNFILTER)]]\
                                           .value_counts().astype(int)
        adreport.fillna(0)
        adreport['CCSfrac']   = (adreport.mappedCCS / adreport.ZMWcount).fillna(0)\
                                                                        .map(PCTFUNC)
        adreport.to_csv('{}/adapterReport.tsv'.format(args.outDir),
                          sep='\t',
                          float_format=FLOATFORMAT)
    #count up read classes
    classCounts = onTarget.groupby(grpCols+enzCols)\
                          .size().rename('ZMWcount')\
                          .reset_index(enzCols)\
                          .fillna(0)
    #add fractions column
    classCounts['Fraction'] = (1.0*classCounts.ZMWcount / classCounts.ZMWcount.sum())\
                              .map(PCTFUNC)

    #write results
    classCounts.to_csv('{}/restrictionCounts.csv'.format(args.outDir),
                         sep='\t',
                         float_format=FLOATFORMAT)

    return classCounts

def adapterParser(adapterFasta, minAlnScore=0,minFlankScore=0):
    '''id,alnScore,flankScore'''
    adMap = {str(i):rec.name for i,rec in enumerate(FastaReader(adapterFasta))}
    def parser(adSeries):
        split = adSeries.str.split(',',expand=True)\
                        .fillna(-1)\
                        .astype({1:float,2:int})
        split.loc[((split[1] < minAlnScore) | (split[2] < minFlankScore)),0] = None
        return split[0].map(adMap)\
                       .fillna(NOADAPTER)
    return parser

def classifyZmw(adapterFasta):
    ads = tuple(sorted(rec.name for rec in FastaReader(adapterFasta)))
    def getClass(data):
        #If any subread has both adapters, call it that
        if np.any(data == ads):
            return ads
        #else if it has subreads with two of the same, call it that
        #possible missing case:  different subreads in same zmw with 
        #ads[0]--ads[0] and ads[1]--ads[1]
        elif np.any(data == (ads[0],ads[0])):
            return (ads[0],ads[0])
        elif np.any(data == (ads[1],ads[1])):
            return (ads[1],ads[1])
        #finally single-sided if nothing else
        else:
            return data[data.map(lambda v:'NoAd' in v)].iloc[0]
    return getClass

def readRestrictionTsv(reFile):
    return pd.read_table(reFile,sep='\s+',names=['enzName','cutsite'])

def getREpos(cutsite):
    site = re.compile(cutsite)
    def rePos(seq):
        return [(m.start(),m.end()) for m in site.finditer(seq)]
    return rePos

class NoAmpRestrictionDiagnostics_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='NoAmpRestrictionDiagnostics.py', description='Generate restriction enzyme and [optionally] adapter report for NoAmp protocol')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='input BAM of CCS alignments to reference')
    parser.add_argument('inBED', metavar='inBED', type=str,
                    help='''input BED defining location of repeat region(s). 
                            Columns ['ctg','start','end','name']''')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta used for mapping BAM')
    parser.add_argument('reTable', metavar='reTable', type=str,
                    help='''tsv table of restriction enzymes and cut sequences.
                            Use multiple lines for degenerate cutsites
                            Columns: ['name','cutsite']''')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='Output directory.  Default %s' % os.getcwd())
    parser.add_argument('-s,--subreadsBAM', dest='subreadsBAM', type=str, default=None,
                    help='''refarmed subreads bam for identifying adapter types. 
                            Must have 'ad' tag from '--adpqc' option in bam2bam.  
                            Default None''')
    parser.add_argument('-a,--adapterFasta', dest='adapterFasta', type=str, default=None,
                    help='fasta of re-farmed adapter sequences.  Required if \'-s\' option used.  Default None')
    parser.add_argument('-m,--minAdapterScore', dest='minAdapterScore', type=float, default=0,
                    help='minimum adapter score.  Default 0')
    parser.add_argument('-f,--minFlankScore', dest='minFlankScore', type=int, default=0,
                    help='minimum adapter flanking score.  Default 0')
    
    try:
        result = main(parser)
    except NoAmpRestrictionDiagnostics_Exception,e:
        print 'ERROR: %s' % e
        sys.exit(1)
