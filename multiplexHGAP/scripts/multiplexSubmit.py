#! python

import json,subprocess,time
import argparse,os,sys,re
from cStringIO import StringIO
from pbcommand.services import ServiceAccessLayer

PRESETS_TEMPLATE='/path/to/presets_template.json'
PBSERVICE       ='/path/to/smrtlink/installdir/smrtcmds/bin/pbservice'
DEFAULTHOST     ='smrtlink'
DEFAULTPORT     =8081
NAMEPOSTFIX     ='_HGAP'
REMOVECHARS     ='();,= '
JOBCSVNAME      ='started_jobs.csv'
    
typeMap = {'float'  :float,
           'integer':int,
           'string' :str,
           'boolean':bool} 

def main(parser,options):
    '''options is a list of json task options added to the parser'''    

    args    = parser.parse_args()
    #must have at least one input
    assert(args.subreadSetID or args.subreadSetIdCsv),'Must define -s or -S'
    #load settings from tempalte
    presets = json.load(open(PRESETS_TEMPLATE))
    
    if args.subreadSetID:
        name  = ''
        ssIdx = {name:int(args.subreadSetID)}
    else:
        ssIdx = parseSubreadsetIdCsv(args.subreadSetIdCsv)

    sal = ServiceAccessLayer(args.host,args.port)

    #prepare file to report jobs started
    columns = ['host','jobId','jobName','jobPath']
    csvfmt  = ','.join(map('{{{}}}'.format,columns))+'\n'
    csvFile = open('{d}/{f}'.format(d=args.outDir,f=JOBCSVNAME),'w')
    #write header
    csvFile.write(','.join(columns)+'\n')
    print 'starting jobs for {i} subreadsets'.format(i=len(ssIdx))
    for name,ssId in ssIdx.items():
        #get the subset
        ss = sal.get_subreadset_by_id(ssId)
        #set the job name
        if args.jobName:
            jobName = args.jobName
        elif name:
            jobName = name + NAMEPOSTFIX
        else:
            jobName = ss['name'] + NAMEPOSTFIX
        presets['name'] = jobName
        #set entry subreadset
        setEntryPoint(presets,ssId)
        #set all options
        for opt in options:
            setTaskOption(presets,opt,getattr(args,opt))         
        #write preset_json
        job_pre = '{d}/{name}_presets.json'.format(d   =os.path.abspath(args.outDir), 
                                                   name=cleanName(jobName.replace(' ','_')))
        with open(job_pre,'w') as oFile:
            json.dump(presets,oFile,indent=2)    
        
        #start job 
        print 'Starting job {name}, {time}'.format(name=jobName,
                                                   time=time.asctime(time.localtime()))
        job = startJob(job_pre,host=args.host,port=args.port)
        jobSummary = job['JOB SUMMARY']
        csvFile.write(csvfmt.format(host   =args.host,
                                    jobId  =int(jobSummary['id']),
                                    jobName=jobSummary['name'],
                                    jobPath=jobSummary['path']))
        if args.wait:
            print 'waiting %i minutes' % args.wait
            time.sleep(60*args.wait)
    csvFile.close()
    
    return None

def parseSubreadsetIdCsv(filename):
    lines   = [line.strip().split(',') 
               for line in open(filename).read().strip().split('\n')]
    jobDict = {name:int(ssId) 
               for name,ssId,jobState in lines
               if ssId.isdigit() and jobState=='SUCCESSFUL'}
    return jobDict

def setEntryPoint(presets,ssId):
    presets['entryPoints'][0]['datasetId'] = ssId
    return presets

def setTaskOption(presets,option,value):
    opt   = filter(lambda p: p['optionId'].endswith(option),presets['taskOptions'])[0]
    dtype = opt['optionTypeId'].split('.')[-1] 
    opt['value'] = typeMap[dtype](value)
    return presets

def cleanName(name,removeChars=REMOVECHARS):
    return re.sub('['+removeChars+']','',name)

def startJob(presetJson,host=DEFAULTHOST,port=DEFAULTPORT):
    cmd  = [PBSERVICE,'run-analysis',presetJson,'--host',host,'--port',port]
    proc = subprocess.Popen( map(str,cmd), stdout=subprocess.PIPE )
    res  = StringIO( proc.communicate()[0] )
    return parseOutput(res)

def parseOutput(stream):
    jobFields = ['DATASET SUMMARY:','JOB SUMMARY:']
    out       = {}
    fld       = None
    for line in stream.read().strip().split('\n'):
        if line in jobFields:
            fld = line[:-1]
            out[fld] = {}
        if line not in jobFields and ':' in line:
            k,v = line.strip().split(':',1)
            out[fld][k.strip()] = v.strip()
    return out

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='multiplexSubmit.py', description='Submit multiple SMRT-Link jobs using json template')
    parser.add_argument('-s,--subreadSetID', dest='subreadSetID', type=int, default=None,
                    help='entry point subreadset')
    parser.add_argument('-S,--subreadSetIdCsv', dest='subreadSetIdCsv', type=str, default=None,
                    help='csv table of entry point subreadsets. columns: {subsetName,subsetId,jobState}')
    parser.add_argument('-o,--outDir', dest='outDir', type=str, default=os.getcwd(),
                    help='locations for job preset files.  default {d}'.format(d=os.getcwd()))
    parser.add_argument('-j,--jobName', dest='jobName', type=str, default=None,
                    help='name for smrtlink jobs (default subreadset name + "%s")' % NAMEPOSTFIX)
    parser.add_argument('-w,--wait', dest='wait', type=int, default=0,
                    help='Minutes to wait between job submissions. Default 0 min')
    parser.add_argument('--host', dest='host', type=str, default=DEFAULTHOST,
                    help='SMRTlink server. Default %s' % DEFAULTHOST)
    parser.add_argument('--port', dest='port', type=int, default=DEFAULTPORT,
                    help='SMRTlink server. Default %i' % DEFAULTPORT)
    #add all taskOptions from template
    template = json.load(open(PRESETS_TEMPLATE))
    options  = []
    for opt in template['taskOptions']:
        task,_,option = opt['optionId'].split('.')
        dtype   = typeMap[opt['optionTypeId'].split('.')[-1]]
        default = opt['value']
        parser.add_argument('--%s'%option, dest=option, type=dtype, default=default,
                            help='{task} {option}. Default {default}'.format(task=task,option=option,
                                                                             default='"%s"'%default if dtype==str
                                                                                     else default))
        options.append(option)

    main(parser,options)
