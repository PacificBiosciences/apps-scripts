#! /usr/bin/python

import sys
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from resources.extract import extractRegion
from resources.utils import countMotifs, \
                            formatMixFloat
from resources.clust import getCounts, \
                            clusterStats, \
                            addHPtag

DEFAULTKMER  = 3
DEFAULTPREFIX= './cluster'
LENGTHFIELD  = 'totalBp'
SLCOLNAMES   = ['Item','Metric','cluster0','cluster1'] 
DEFAULTCI    = [0.025,0.975] #95%
AGGFUNCS     = [np.median,np.mean]
FLOATFMT     = '%.1f'

def main(parser):
    args = parser.parse_args()
    
    motifs       = args.motifs.split(',')
    motifCounter = countMotifs(motifs,
                               lengthField=LENGTHFIELD,
                               collapseHP=args.collapseHP)    
    seqGen       = extractRegion(args.inBAM,
                                 args.reference,
                                 region=args.region,
                                 flanksize=args.flanksize,
                                 revcomp=args.revcomp)

    print("Reading sequence")
    kmerCounts,motifCounts = getCounts(seqGen,args.kmer,motifCounter) 
    if not len(kmerCounts) > 1:
        print(f'No Records in {args.inBAM}. Must have > 1 records to cluster')
        return None
    
    #kmeans clustering
    print("Clustering")
    kmeans = KMeans(n_clusters=args.clusters)
    kmeans.fit(kmerCounts)
    clusterIdx = kmeans.predict(kmerCounts)

    #get cluster stats
    columns       = motifs + [LENGTHFIELD]
    names,results = clusterStats(motifCounts,clusterIdx,columns,
                                 aggFuncs=AGGFUNCS,randomSeed=args.seed,
                                 ci=DEFAULTCI)

    #write results
    print("Writing Results")
    #results.to_excel('%s.summary.xlsx' % args.prefix)
    results.to_csv(f'{args.prefix}.summary.csv', float_format=FLOATFMT)

    if args.smrtlink:
        #add blank column if nclusters set to 1 by user
        if args.clusters==1:
            results.loc['cluster1'] = None
        #Transpose results table and set colnames
        resultsSL = results.T.reset_index().set_axis(SLCOLNAMES,axis=1)
        #Remove duplicate "Items"/first column for esthetics
        itCol = SLCOLNAMES[0]
        resultsSL[itCol][resultsSL[itCol].duplicated(keep='first')] = ''
        #Better counts field name for SL
        resultsSL[itCol][resultsSL[itCol] == 'Read'] = 'Hifi Reads'
        #format floats in mixed columns
        for col in SLCOLNAMES[-2:]:
            resultsSL[col] = resultsSL[col].map(formatMixFloat(FLOATFMT))
        #write it
        resultsSL.to_csv(f'{args.prefix}.summarySL.csv',
                         index=False)

    with open(f'{args.prefix}.readnames.txt', 'w') as namefile:
        for cluster,reads in motifCounts.groupby(names.reindex(clusterIdx).values):
            namefile.write(f'>{cluster}\n')
            namefile.write('\n'.join(reads.index) + '\n')

    if not args.noBam:
        #strip extra fields from names
        readNames  = motifCounts.index.str.split('/').str[:3].str.join('/')
        clusterMap = dict(list(zip(readNames,clusterIdx)))
        outBam     = f'{args.prefix}.hptagged.bam'
        addHPtag(args.inBAM,outBam,clusterMap,dropNoClust=args.drop)

    print("Done")

    return kmerCounts,motifCounts,results

class ClusterByRegion_Exception(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='clusterByRegion.py', description='kmer clustering by target region')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='input BAM of CCS alignments')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta used for mapping BAM.  Must have .fai index.')
    parser.add_argument('region', metavar='region', type=str,
                    help='Target region, format \'[chr]:[start]-[stop]\'.  Example \'4:3076604-3076660\'')
    parser.add_argument('-m','--motifs', dest='motifs', type=str, default=None,required=True,
                    help='comma-separated list of motifs to count')
    parser.add_argument('-k','--kmer', dest='kmer', type=int, default=DEFAULTKMER,
                    help=f'kmer size for clustering. Default {DEFAULTKMER}')
    parser.add_argument('-c','--clusters', dest='clusters', type=int, default=2,
                    help='clusters/ploidy count. Default 2')
    parser.add_argument('-r','--revcomp', dest='revcomp', action='store_true',
                    help='Reverse complement extracted sequence')
    parser.add_argument('-p','--prefix', dest='prefix', type=str, default=DEFAULTPREFIX,
                    help=f'Output prefix. Default {DEFAULTPREFIX}')
    parser.add_argument('-f','--flanksize', dest='flanksize', type=int, default=100,
                    help='Size of flanking sequence mapped for extracting repeat region.  Default 100')
    parser.add_argument('-s','--seed', dest='seed', type=int, default=42,
                    help='Seed for resampling ci95.  Default 42')
    parser.add_argument('-x','--noBam', dest='noBam', action='store_true',
                    help='Do not export HP-tagged bam of clustered reads')
    parser.add_argument('--smrtlink', dest='smrtlink', action='store_true',
                    help='Export smrtlink-formatted report')
    parser.add_argument('-d','--drop', dest='drop', action='store_true',
                    help='Drop reads with no cluster in output bam.  Default keep all reads.')
    parser.add_argument('-u','--collapseHP', dest='collapseHP', action='store_true',
                    help='Collapse homopolymers before analysis.  Default use original sequence.')

    try:
        main(parser)
    except ClusterByRegion_Exception as e:
        print(f'ERROR: {e}')
        sys.exit(1)
