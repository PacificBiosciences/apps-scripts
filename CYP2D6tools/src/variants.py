#!/usr/bin/env python3
__version__ = '0.1.0'
import re,pysam,os
import mappy as mp
import pandas as pd
import hashlib
from operator import itemgetter
from collections import Counter
from . import config as cfg

class Caller:
    def __init__(self,consensusFastas,runName,reference,spacer,alnPreset,sMap,minFrac=0.01,
                 ignoreMissing=False,read_info=None,readsFq=None,datetime=None):
        self.inputFa                = consensusFastas
        self.runName                = runName
        self.aligner                = Aligner(reference,preset=alnPreset)
        self.spacerAln              = Aligner(spacer,preset=alnPreset)
        self.sMap                   = sampleMap(sMap)
        self.minFrac                = minFrac
        self.ignoreMissing          = ignoreMissing
        self.datetime               = datetime
        self.alleles, self.variants = self.callVariants()
        #self.
        if read_info and readsFq:
            self.getSupport(read_info,readsFq)

    def callVariants(self):
        alleles,variants = [],[]
        for consensusFa in self.inputFa:
            consensusType = self._parsePbAAfastaName(consensusFa)
            for rec in pysam.FastxFile(consensusFa):
                aln = self.aligner(rec,skipFailed=(consensusType == 'failed'))
                if aln is None:
                    continue #skip failed consensus if they do not map
                nameDict  = self._parseRecordName(rec)
                if consensusType == 'failed' and nameDict.get('cluster_freq',0) < self.minFrac:
                    continue #skip it
                bioSample = self.sMap(nameDict['barcode'])
                tableKey  = self._getKey(rec.name,bioSample,self.runName,os.path.abspath(consensusFa))
                nameDict.update({'uuid'         :tableKey,
                                 'bioSample'    :bioSample,
                                 'runName'      :self.runName,
                                 'source'       :os.path.abspath(consensusFa),
                                 'clusterStatus':consensusType,
                                 'chrom'        :aln.ctg,
                                 'alnStart'     :aln.r_st,
                                 'alnStop'      :aln.r_en,
                                 'hasSpacer'    :self.hasSpacer(rec),
                                 'csString'     :aln.cs,
                                 'length'       :len(rec.sequence),
                                 'datetime'     :self.datetime})
                alleles.append(pd.Series(nameDict))
                variants.append(self.makeVarTable(aln,tableKey))

        #check for missing/failed results
        if len(self.sMap.remaining) and not self.ignoreMissing:
            for bc,sample in self.sMap.remaining:
                print(f'WARNINING: Failed sample {sample} / {bc}')
                tableKey = self._getKey(bc,sample)
                alleles.append(pd.Series({'uuid'     :tableKey,
                                          'bioSample':sample,
                                          'barcode'  :bc}))
        if len(alleles) == 0: #no recs
            alleles_out = self._getEmpty(cfg.database['alleleTable']).set_index('uuid')
            variant_out = self._getEmpty(cfg.database['variantTable']).set_index(['uuid','CHR','POS'])
        else: 
            alleles_out = pd.DataFrame(alleles).set_index('uuid')
            variant_out = pd.concat(variants)
        return alleles_out,variant_out 

    def makeVarTable(self,aln,key):
        variants = pd.DataFrame(parseCS(aln.cs,aln.r_st),columns=['POS','VAR'])
        variants['CHR'] = aln.ctg
        variants['uuid'] = key
        variants.set_index(['uuid','CHR','POS'],inplace=True)
        return variants

    def hasSpacer(self,rec):
        aln = self.spacerAln(rec)
        #TODO add more sanity checks here
        return int(aln is not None)

    def addRefCalls(self,variants,alleles):
        #query to find ref-call consensus reads that span variants 
        qry = 'chrom==@ctg \
             & alnStart < @pos <= alnStop \
             & uuid not in @vnts.index.get_level_values("uuid")'
        try:
            refCalls = pd.DataFrame([{'uuid': uuid,
                                      'CHR' : ctg,
                                      'POS' : pos,
                                      'VAR' : '.'}
                                     for (ctg,pos),vnts in variants.groupby(['CHR','POS'])
                                     for uuid in alleles.query(qry).index])\
                         .set_index((['uuid','CHR','POS']))    
            return pd.concat([variants,refCalls])
        except KeyError: #no records to add
            return variants
    
    def getSupport(self,read_info,readsFq):
        readSupport   = supportEngine(read_info,readsFq,self)
        vnts          = self.addRefCalls(self.variants,self.alleles)
        self.variants = readSupport(self.alleles,vnts)
        return None

    def _parseRecordName(self,rec):
        resDict = re.compile(cfg.caller['namePattern']).search(rec.name).groupdict()
        for k in ['numreads','cluster']:
            resDict[k] = int(resDict[k])
        splits  = list(map(lambda s:s.strip().split(),rec.comment.split(':')))
        flds    = map(itemgetter(-1),splits[:-1])
        vals    = [self._safeFloat(s[0]) for s in splits[1:-1]] + [' '.join(splits[-1])]
        resDict.update(dict(zip(flds,vals)))
        return resDict

    def _parsePbAAfastaName(self,fa):
        if fa.endswith('passed_cluster_sequences.fasta'):
            return 'passed'
        if fa.endswith('failed_cluster_sequences.fasta'):
            return 'failed'
        else:
            raise Variants_Error(f'Input fasta {fa} not in pbAA format')

    def _safeFloat(self,v):
        try:
            return float(v)
        except ValueError:
            return v

    def _getKey(self,*args):
        return hashlib.md5(''.join(map(str,args)).encode()).hexdigest()

    def _getEmpty(self,tbl):
        return pd.DataFrame(columns=cfg.tableMap[tbl].keys())

class Aligner:
    presets = {'splice'    : {'preset' :'splice'},
               'map-pb'    : {'preset' :'map-pb'},
               'gaplenient': {'scoring':(1,2,2,1,18,0)},  # (A,B,o,e,O,E)
               'gapstrict' : {'scoring':(2,5,10,4,56,1)}} # equiv to pbmm2 align --preset CCS -o 10

    def __init__(self,reference,preset='gaplenient'):
        self.kwargs = {'fn_idx_in' : reference,
                       'best_n'    : 1}
        self.kwargs.update(self.presets[preset])
        self._aligner = mp.Aligner(**self.kwargs)
    
    def __call__(self,rec,skipFailed=True):
        try:
            return list(filter(lambda a:a.is_primary,self._aligner.map(seq=rec.sequence,cs=True)))[0]
        except IndexError:
            if skipFailed:
                return None
            else:
                raise pbCYP2D6_Error(f'Unable to align {rec.name}')

def parseCS(csString,start=0,zeroIndex=True):
    ops = ':*-+~' 
    op  = None
    val = ''
    i   = start + int(zeroIndex)
    for s in csString:
        if s in ops:
            if op:
                incr,vnt = parseOp(i,op,val)
                if op=='+': #ins, place on prev position
                    i -= 1
                if vnt:
                    yield i,vnt
                i += incr
                val = ''
            op = s
        else:
            val += s
    incr,vnt = parseOp(i,op,val)
    if vnt:
        yield i,vnt

_splice = re.compile('[acgtn]{2}([0-9]+)[acgtn]{2}')
def parseOp(i,op,val):
    if op == ':': #match
        return int(val),None
    elif op == '-': #del
        return len(val),op+val
    elif op == '+': #ins
        return 1,op+val
    elif op == '*': #mismatch
        return 1,op+val
    elif op == '~': #splice/large del
        try:
            size = int(_splice.search(val).groups()[0])
        except Exception:
            raise Variants_Error(f'Failed to parse variant: {op+val}')
        return size,f'{op}{size}del'
    else:
        raise Variants_Error(f'Unknown variant type: {op}')

class supportEngine:
    FIELDNAME=cfg.database['supportField']

    def __init__(self,readInfo,hifiReads,vcaller):
        self.readInfo  = self._loadInfo(readInfo)
        self.vcaller   = vcaller
        self.alnGen    = self._align(hifiReads)
        self.vTable    = self._makeTable()
        #for now, only single sample runs
        assert self.vTable.Sample.nunique()==1, "can only do single samples"

    def __call__(self,alleles,variants):
        if variants.empty: #quick skip to avoid key error
            return variants.assign(**{self.FIELDNAME:None})

        sVar = pd.merge(alleles[['guide','cluster','numreads']],
                        variants,
                        left_index=True,
                        right_index=True)
        sVar.cluster = sVar.cluster.astype(int)

        bigTable = pd.merge(sVar.reset_index(),
                            self.vTable[['CHR','POS','VAR','guide','ClusterId']],
                            how='left',
                            left_on=['CHR','POS','guide','cluster'],
                            right_on=['CHR','POS','guide','ClusterId'],
                            suffixes=['_call','_reads'])

        counts   = bigTable.groupby(['uuid','CHR','POS'])\
                           .apply(self._countVars)\
                           .rename(self.FIELDNAME)

        return variants.join(counts)

    def _loadInfo(self,readInfo):
        return pd.read_csv(readInfo,sep=' ',names=cfg.caller['readInfo'])\
                 .set_index('readName')[['Sample','guide','ClusterId']]

    def _align(self,hifi):
        for rec in pysam.FastxFile(hifi):
            aln = self.vcaller.aligner(rec)
            if aln is not None:
                yield rec.name,aln

    def _countVars(self,reads):
        #numreads = int(reads.numreads.iloc[0])
        numreads = len(self.vTable.query(f'guide==@reads.guide.iloc[0] \
                                         & ClusterId==@reads.ClusterId.iloc[0] \
                                         & alnStart < @reads.POS.iloc[0] < alnStop').index.unique())
        if reads.VAR_reads.isnull().all():
            return {'.':numreads}
        cnts     = Counter(reads.VAR_reads)
        subtot   = sum(cnts.values())
        if subtot < numreads:
            cnts['.'] += numreads - subtot
        return dict(cnts)

    def _makeTable(self):
        vTable = pd.concat([self.vcaller.makeVarTable(aln,readName)\
                             .assign(alnStart=int(aln.r_st),alnStop=int(aln.r_en))
                            for readName,aln in self.alnGen])\
                   .reset_index(['CHR','POS'])

        #identify deletions >1 base
        #insert '*' if deleted positions are variants elsewhere
        #to prevent the reads from being counted as refcall over those pos
        multiDel     = vTable[vTable.VAR.str.startswith('-') & (vTable.VAR.str.len() > 2)]
        varPositions = set(map(tuple,vTable[['CHR','POS']].drop_duplicates().values))
        try:
            deletions    = pd.DataFrame([{'readName':row.name,
                                          'CHR'     :row.CHR,
                                          'POS'     :pos,
                                          'VAR'     :'*'}
                                          for i,row in multiDel.iterrows()
                                          for pos in range(row.POS+1,row.POS+len(row.VAR)-1)
                                           if (row.CHR,pos) in varPositions ])\
                                         .set_index('readName')
            return pd.concat([vTable,deletions]).join(self.readInfo)
        except KeyError:  #nothing to add
            return vTable.join(self.readInfo)

class sampleMap:
    _UNK='unknown_sample'

    def __init__(self,sampleMapFile=None):
        self.sMap      = {}
        self.remaining = []
        self.mapfile   = sampleMapFile
        #self.__call__  = self._getCallFunc(sampleMapFile)
        self.callFunc  = self._getCallFunc(sampleMapFile)

    def __call__(self,bc):
        return self.callFunc(bc)

    def _getCallFunc(self,mapFile):
        if mapFile is None:
            return (lambda bc: 'None')
        else:
            self.sMap      = dict(pd.read_csv(mapFile)[cfg.caller['sampleMapCols']].values)
            self.remaining = list(self.sMap.items()) 
            return self._getSample

    def _getSample(self,bc):
        #sn = self.sMap[bc]
        sn = self.sMap.get(bc,self._UNK)
        if sn == self._UNK:
            print(f'WARNING:  Barcode {bc} not in sample map {self.mapfile}')
        try:
            self.remaining.pop(self.remaining.index((bc,sn)))
        except ValueError:
            pass
        return sn

class Variants_Error(Exception):
    pass
