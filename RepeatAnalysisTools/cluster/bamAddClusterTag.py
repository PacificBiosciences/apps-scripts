#! /usr/bin/python

from resources.clust import readClusterFile,addHPtag

def main(parser):
    args = parser.parse_args()    

    readMap = readClusterFile(args.clusters)
    addHPtag(args.inBam,args.outBAM,readMap,dropNoClust=args.drop)

    return None

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='bamAddClusterTag.py', description='Update BAM so that HP tag matches cluster')
    parser.add_argument('inBAM', metavar='inBAM', type=str,
                    help='BAM file to edit')
    parser.add_argument('clusters', metavar='clusters', type=str,
                    help='Cluster file')
    parser.add_argument('-o,--outBAM', dest='outBAM', type=str, default=sys.stdout,
                    help='Output BAM file.  Default stdout' )
    parser.add_argument('-d,--drop', dest='drop', action='store_true',
                    help='Drop record with no defined cluster.  Default write all records' )

    main(parser)
