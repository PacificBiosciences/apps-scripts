#! /home/UNIXHOME/jharting/anaconda2/bin/python

import os,sys,argparse
import pandas as pd
from pbcommand.services import ServiceAccessLayer

DEFAULTHOST= 'smrtlink-alpha'
DEFAULTPORT= 8081
DEFAULTCSV = 'collated_reports'
FINISHED   = ['SUCCESSFUL','FAILED','TERMINATED']
LINKFMT    = 'http://{0[0]}:8080/#/analysis/job/{0[1]}'.format
#FLOATFMT   = '{:.1f}'.format


def main(parser):

    args = parser.parse_args()

    jobs = pd.read_csv(args.jobCsv)
    sal  = ServiceAccessLayer(args.host,args.port)
    #get dicts of values for all jobs
    rpts = jobs.jobId.apply(sal.get_analysis_job_report_attrs).values
    #check for unfinished jobs, exit if any
    unfinished = [j for j in jobs.jobId
                  if sal.get_job_by_id(j).state not in FINISHED] 
    if unfinished:
        for j in unfinished:
            print 'job %i still running'%j
        print 'Exiting'
        sys.exit()
    #put the reports together and index with (jobName,host,jobId,link)
    jobs['link'] = jobs[['host','jobId']].apply(LINKFMT,axis=1) 
    columns      = ['jobName','host','jobId','link']
    idx          = pd.MultiIndex.from_arrays(jobs[columns].values.T,names=columns)
    collated     = pd.DataFrame.from_records(rpts,index=idx).T
    for fmt,fnc in zip(['.csv','.xls'],
                       [pd.DataFrame.to_csv,pd.DataFrame.to_excel]):
        ofile = '{d}/{name}{fmt}'.format(d=args.outDir,
                                         name=DEFAULTCSV,
                                         fmt=fmt)
        fnc(collated,ofile)  # float_format=FLOATFMT
        print 'Wrote results to %s' % ofile
    
    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='collateReports.py', description='put together reports listed in a csv by jobId,jobName')
    parser.add_argument('jobCsv', metavar='jobCsv', type=str,
                    help='csv file listing jobs.  columns must include: "jobId,jobName"')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='directory for output. Default %s' % os.getcwd())
    parser.add_argument('--host', dest='host', type=str, default=DEFAULTHOST,
                    help='SMRTlink server.  Default %s' % DEFAULTHOST)
    parser.add_argument('--port', dest='port', type=int, default=DEFAULTPORT,
                    help='SMRTlink server. Default %i' % DEFAULTPORT)
    

    main(parser)
