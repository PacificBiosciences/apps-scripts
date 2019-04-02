import pysam
import sys

def makeID(*args):
    return '_'.join(map(str,args))

def main(parser):
    args = parser.parse_args()

    if args.barcodes:
        bcMap = {i:rec.name for i,rec in enumerate(pysam.FastxFile(args.barcodes))}
    else:
        sys.stderr.write('INFO: barcodes=None, using default bc1001 ... bc1384\n')
        bcMap = {i:'bc1{:03}'.format(i+1) for i in range(384)}

    makeBarcode = lambda bcTag: '{}--{}'.format(*[bcMap[i] for i in bcTag])

    with pysam.AlignmentFile(args.infile) as infile:
        #Copy header from infile
        header = infile.header.as_dict()
        #get existing RG tags as templates
        existing = {ds['ID']:ds for ds in header['RG']}
        #reset header RG elements
        header['RG'] = []
        #loop through infile to get all unique barcode/RG combos
        uniqueTags = {(makeBarcode(rec.get_tag('bc')),rec.get_tag('RG')) for rec in infile}
        #create an RG row for each pair and add to header
        for barcode,rgId in uniqueTags:
            new = existing[rgId].copy()
            new['ID'] = makeID(barcode,rgId)
            new['SM'] = barcode
            header['RG'].append(new)
        #write rows from infile after editing RG tag
        with pysam.AlignmentFile(args.outfile,'wb',header=header) as outfile:
            infile.reset()
            for rec in infile:
                bc  = makeBarcode(rec.get_tag('bc'))
                idx = rec.get_tag('RG')
                rec.set_tag('RG',makeID(bc,idx))
                outfile.write(rec)
    return None


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='updateBamBarcode2SM.py', description='Update RG tags with barcode information, and setting SM=barcode')
    parser.add_argument('infile', metavar='infile', type=str,
                    help='BAM file to edit')
    parser.add_argument('-b,--barcodes', dest='barcodes', type=str, default=None,
                    help='No header in m5 file.  Default use standard PB barcodes, lexicographically ordered: bc1001..bc1384' )
    parser.add_argument('-o,--outfile', dest='outfile', type=str, default=sys.stdout,
                    help='Output BAM file.  Default stdout' )

    main(parser)
