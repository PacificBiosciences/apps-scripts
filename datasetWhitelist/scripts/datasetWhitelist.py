from pbcore.io import DataSet,Filters,FastaReader

def main(parser):
    args = parser.parse_args()

    filt  = Filters()
    dset  = DataSet(args.inXml)
    reads = FastaReader(args.inFasta)
    if args.subreads:
        filt.addRequirement(QNAME=[('=',rec.name) for rec in reads])
    else:
        filt.addRequirement(zm=[('=',hn) for hn in set(map(getZmw,reads))])
    dset.addFilters(filt)
    dset.write(args.outXml)

def getZmw(rec):
    return int(rec.name.split('/')[1])

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='datasetWhitelistZmw.py', description='Generate a whitelisted dataset xml from an input fasta of reads (subreads or ccs)')
    parser.add_argument('inXml', metavar='inXml', type=str,
                    help='input pacbio subread dataset.')
    parser.add_argument('inFasta', metavar='inFasta', type=str,
                    help='fasta of subread or ccs sequences.')
    parser.add_argument('outXml', metavar='outXml', type=str,
                    help='output xml.' )
    parser.add_argument('-s,--subreads', dest='subreads', action='store_true', default=False,
                    help='whitelist subread names instead of zmws (will only work with fasta of subreads).  default false.')

    main(parser)

