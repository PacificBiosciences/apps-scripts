#! python

from pbcore.io import SubreadSet
from pbcommand.services import ServiceAccessLayer
import sys,os,argparse,time,re

BCSEP      ='--'
NAMEFMT    ='{prefix} (barcode={barcode})'
#NAMEFMT='{prefix} ({barcode};filter={filter})'
PATHFMT    ='{dname}/{prefix}_{barcode}.subreadset.xml'
REMOVECHARS='();,= '
FINISHED   = ['SUCCESSFUL','FAILED','TERMINATED']

DEFAULTHOST='smrtlink'
DEFAULTPORT=8081
DEFAULTCSV ='uploaded_subreadsets.csv'

def main(parser):
    
    args = parser.parse_args()
    #connect to service
    sal  = ServiceAccessLayer(args.host,args.port)
    #if barcode fasta, use those names, else the index
    if args.barcodeFasta:
        from pbcore.io import FastaReader
        barcodes = [rec.name for rec in FastaReader(args.barcodeFasta)]
        f = lambda i: barcodes[i]
    else:
        f = str
    def getBcName(idxPair):
        return BCSEP.join(map(f,idxPair))
    
    sset = SubreadSet(args.subreadSet)
    if args.namePrefix:
        prefix = args.namePrefix
    else:
        prefix = sset.name
    #get rid of any spaces in name
    prefix = cleanName(prefix.replace(' ','_'))

    #columns in csv
    outfmt = '{subsetName},{subsetId},{jobState}\n'
    #out csv
    oFile  = open('{d}/{f}'.format(d=args.outDir,f=args.outCsvName),'w')
    #write header
    oFile.write(outfmt.format(subsetName='subsetName',subsetId='subsetId',jobState='jobState'))
    print "Splitting by barcode"
    bcSets  = sset.split(barcodes=True)
    print "Importing {n} barcodes to {host}:{port}".format(n=len(filter(lambda s:len(s)>=args.minSubreads,bcSets)),
                                                           host=args.host,
                                                           port=args.port)
    for bcSet in bcSets:
        barcode    = getBcName(eval(bcSet.barcodes[0]))
        #check for minimum number of subreads or skip
        if bcSet.numRecords < args.minSubreads:
            print 'skipping barcode {barcode}. Too few subreads ({numRecs})'.format(barcode=barcode,
                                                                                    numRecs=bcSet.numRecords)
            continue
        #add bc quality filter
        bcSet.filters.addRequirement(bq=[('>=',args.filter)])        
        #set name (assume one and only one barcode)
        bcSet.name = NAMEFMT.format(prefix=prefix,barcode=barcode,filter=args.filter)
        #xml name
        bcSetPath  = PATHFMT.format(dname=os.path.abspath(args.outDir),prefix=prefix,barcode=barcode)
        bcSet.write(bcSetPath)
        #import -- check if already imported first
        ss = sal.get_dataset_by_uuid(bcSet.uuid)
        if ss:
            state = 'SUCCESSFUL'
            print 'dataset {name} previously imported'.format(name=bcSet.name)
        else:
            job = sal.import_dataset_subread(bcSetPath)
            print 'waiting for import job {i}...'.format(i=job.id)
            while job.state not in FINISHED:
                time.sleep(3)
                job = sal.get_job_by_id(job.id)
            state = job.state
            print 'job {i}: import {name}: {state}'.format(i    =job.id,
                                                           name =bcSet.name,
                                                           state=state)
        if state == 'SUCCESSFUL':
            ss = sal.get_dataset_by_uuid(bcSet.uuid)
            oFile.write(outfmt.format(subsetName=ss['name'],subsetId=ss['id'],jobState=state))
    oFile.close()
    return None

def cleanName(name,removeChars=REMOVECHARS):
    return re.sub('['+removeChars+']','',name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='splitBarcodeUpload.py', description='split barcoded subreadset and import subsets to server')
    parser.add_argument('subreadSet', metavar='subreadSet', type=str,
                    help='barcoded dataset xml')
    parser.add_argument('-b,--barcodeFasta', dest='barcodeFasta', type=str, default=None,
                    help='fasta of barcode sequences used for barcoding.  Will be used to name outputs if given')
    parser.add_argument('-f,--filter', dest='filter', type=int, default=45,
                    help='Min bcQual filter. Default 45' )
    parser.add_argument('-r,--minSubreads', dest='minSubreads', type=int, default=5000,
                    help='Min number of subreads per barcode. Default 5000' )
    parser.add_argument('-n,--namePrefix', dest='namePrefix', type=str, default=None,
                    help='prefix name for uploading (default subreadset name).  Will postfix barcode: [namePrefix]_[barcode]')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='directory to save xml for uploading. Default %s' % os.getcwd())
    parser.add_argument('--outCsvName', dest='outCsvName', type=str, default=DEFAULTCSV,
                    help='csv file of uploaded subreadset ids. Default [outDir]/%s' % DEFAULTCSV)
    parser.add_argument('--host', dest='host', type=str, default=DEFAULTHOST,
                    help='SMRTlink server.  Default %s' % DEFAULTHOST)
    parser.add_argument('--port', dest='port', type=int, default=DEFAULTPORT,
                    help='SMRTlink server. Default %i' % DEFAULTPORT)
    

    main(parser)
