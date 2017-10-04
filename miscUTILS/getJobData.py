#! /usr/bin/python

import os
from shutil import copy
from pbcommand.services import ServiceAccessLayer

DEFAULTHOST='smrtlink'
DEFAULTPORT=8081
#output field size for printing names
fmt        ='{:36}'.format
#data file attributes to print with -n option
#attrs      = ['name','uuid','path']
attrs      = ['name','path']

def main(parser):
    args = parser.parse_args()

    if not args.names:
        assert args.dataName, 'Must pass dataName if not using -n option'

    #set up connection to service
    sal    = ServiceAccessLayer(args.host,args.port)
    #get datastore 
    dstore = sal.get_analysis_job_datastore(args.jobNumber)
    #loop through data
    for uuid,dsfile in dstore.files.items():
        if args.names:
            #print the attribute values
            print '\t'.join([fmt(getattr(dsfile,a))
                             for a in attrs])    
        elif dsfile.name in args.dataName:
            #cp file to outdir
            ofile = '{o}{s}{n}'.format(o=args.outDir,
                                       s=os.path.sep,
                                       n=os.path.basename(dsfile.path))
            copy(dsfile.path,ofile)
            print '\t'.join([fmt('\'%s\''%dsfile.name),'=>',ofile])

    return 

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(prog='getJobData.py', description='Get job outputs from smrtlink datastore')
    parser.add_argument('jobNumber', metavar='jobNumber', type=int,
                    help='SMRT Link job number')
    parser.add_argument('dataName', metavar='dataName', type=str, default=None, nargs='*',
                    help='name of datasetset to get.  default None.  (use -n to see available names)')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='output destination directory. default cwd' )
    parser.add_argument('-n,--names', dest='names', action='store_true', default=False,
                    help='print names of available items in datastore. default False' )
    parser.add_argument('--host', dest='host', type=str, default=DEFAULTHOST,
                    help='smrtlink host name. default %s' % DEFAULTHOST )
    parser.add_argument('--port', dest='port', type=int, default=DEFAULTPORT,
                    help='smrtlink services port. default %i' % DEFAULTPORT )

    main(parser)
