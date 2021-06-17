#!/usr/bin/env python3
__version__ = '0.1.1'
import pysam,sys
import numpy as np
import pandas as pd
from operator import itemgetter
from collections import Counter
try:
    import src.config as cfg
except ModuleNotFoundError:
    import config as cfg


DEFAULTMINFREQ=0.01
DEFAULTQUAL=200
SUPPORTFIELD='supportReads'

def main(parser):
    args = parser.parse_args()

    print('Loading data')
    newVcf = VcfCreator(args.out or '-',
                        args.alleles,
                        args.variants,
                        args.reference,
                        'wb' if args.compress else 'w',
                        args.sampleCol,
                        args.minFreq,
                        args.passOnly,
                        args.merge,
                        args.database,
                        args.query)

    print('Generating Vcf')
    newVcf.run()         
    
    return newVcf


class VcfCreator:
    def __init__(self,filename,alleles,
                 variants,referenceFa,
                 mode='w',sampleCol='barcode',
                 minFreq=DEFAULTMINFREQ,passOnly=False,
                 mergeVars=False,database=None,
                 query=None,dataframe=False):
        self.sampleCol   = sampleCol
        self.mergeVars   = mergeVars
        self.isDataframe = dataframe
        if passOnly:
            self.alleles = self.openAlleles(alleles,database,query).query('clusterStatus == "passed"')
        else:
            self.alleles = self.openAlleles(alleles,database,query).query('cluster_freq >= @minFreq')
        self.variants  = self.openVariants(variants,database,query)
        self.samples   = self.alleles[sampleCol].unique()
        self.filtInput = self.makeFilteredInput()
        self.filename  = filename
        self.mode      = mode
        self.reference = pysam.FastaFile(referenceFa)
        self.header    = self.makeHeader()
        self.vcfFile   = pysam.VariantFile(filename,mode=mode,header=self.header)
        
    def run(self):
        self._loadSamples()
        self.vcfFile.close()
        #if self.filename and self.mode=='wb':
        #    pysam.index(self.filename)
        return self

    def makeFilteredInput(self):
        support = cfg.database['supportField']
        filtInput = self.alleles.join(self.variants).reset_index('POS')
        #check for empty variant set
        if filtInput.empty:
            return filtInput.set_index('POS',append=True)
        #adjust deletion position to include prev base in ref
        #also duplicate ref call for shifted dels
        delPos = filtInput.VAR.str.startswith('-')
        shiftRef = filtInput[filtInput.POS.isin(filtInput[delPos].POS.unique()) & (filtInput.VAR=='.')].copy()
        filtInput.POS -= delPos.astype(int)
        if len(shiftRef):
            shiftRef.POS -= 1
            filtInput = pd.concat([filtInput,shiftRef])        
        if support in filtInput.columns:
            filtInput['nsupport'] = filtInput.apply(lambda r:r[support].get(r['VAR'],0),axis=1)
        else:
            filtInput['nsupport'] = 0
            filtInput[support]    = [{'dummy':0}]*len(filtInput)
        return filtInput.set_index('POS',append=True)

    def openAlleles(self,alleles,database=None,query=None):
        dtype = {'numreads':int}
        if self.isDataframe:
            return alleles
        if database:
            df = self.dbImport(database,alleles,query).set_index('uuid')
            if df.empty:
                raise Pbaa2Vcf_Error(f'No records returned for db "{database}" and query "{query}') 
            for col,t in dtype.items():
                df[col] = df[col].astype(t)
            return df
        else:
            return pd.read_csv(alleles,index_col=0,dtype=dtype)

    def openVariants(self,variants,database=None,query=None):
        support = cfg.database['supportField']
        if self.isDataframe:
            return variants
        if database:
            df = self.dbImport(database,variants,query)      
        else:
            dtype = {'CHR':str,
                     'POS':int}
            df = pd.read_csv(variants,dtype=dtype)
        if support in df.columns:
            df[support] = df[support].map(eval)
        return df.set_index(['uuid','CHR','POS'])
        
    def dbImport(self,database,table,query):
        import sqlalchemy as sqla
        import sqlite3 
        dbengine = sqla.create_engine(f'sqlite:///{database}', echo=False)
        sql = f'SELECT * FROM {table}'
        if query:
            sql = f'{sql} WHERE uuid IN ({query})'
        try:
            return pd.read_sql(sql,con=dbengine)
        except (sqla.exc.OperationalError,sqlite3.OperationalError) as e:
            #raise Pbaa2Vcf_Error(f'Unable to access {table} in {database}')
            raise Pbaa2Vcf_Error(f'Unable to execute query in {database}. DB error as follows:\n\n{e}')

    def makeHeader(self):
        header = pysam.VariantHeader()

        header.filters.add('PBAAFAIL',None,None,'Consensus failed pbAA filters')
        header.filters.add('HP',None,None,'Homopolymer length variant')
        header.info.add('NS',1,'Integer','Number of samples with data')
        header.info.add('AF','A','Float','Allele frequency')
        header.formats.add('GT',1,'String','Genotype ')
        header.formats.add('DP',1,'Integer',"Read depth")
        #maybe
        header.formats.add('FT','.','String','pbAA filter')
        #header.formats.add('GQ',1,'Integer',"Conditional genotype quality")
        header.formats.add('AQ','.','Float',"pbAA mean cluster quality")
        #header.formats.add('MIN_DP',1,'Integer',"Minimum DP observed within the GVCF block.")
        header.formats.add('AD','.','Integer',"Reads supporting each alt call")
        header.formats.add('VAF','.','Float',"pbAA cluster frequency")
        #header.formats.add('PL','G','Integer',"Phred-scaled genotype likelihoods rounded to the closest integer")
        header.formats.add('TG','.','String','pbAA guide')
        header.formats.add('HP','.','Integer',"pbAA cluster identifier")
        header.formats.add('DV','.','Float',"pbAA diversity score")
        header.formats.add('CH','.','Float',"pbAA chimera score")

        for ctg,length in zip(self.reference.references,self.reference.lengths):
            header.contigs.add(ctg,length=length)

        header.add_meta('commandline',value=self._getCommandLine())

        #add samples
        for s in self.samples:
            header.add_sample(s)
        
        return header
            
    def _loadSamples(self):
        for loc,data in self.filtInput.groupby(['CHR','POS']): 
            if (data.VAR == '.').all(): #variants filtered out, but refcall remains
                continue
            rec = self.new_record(loc,data)
            self.vcfFile.write(rec)
            
    def _getCommandLine(self):
        #TODO
        return 'pbaa2vcf <opts> <args>'

    def _getAlleles(self,alts):
        row     = alts.loc[alts.ref.str.len().idxmax()]
        alleles = [row.ref] + sorted(alts.alt.unique()) 
        return alleles

    def _makeAlts(self,data,altFunc):
        alts  = pd.DataFrame.from_records(data.VAR.map(altFunc),
                                          index=data.index).dropna()
        #check for multi var at same locus within sample
        if not alts.index.is_unique:
            alts = alts.groupby(['uuid','CHR','POS']).apply(self._mergeAlts)
        return alts
            
    def _mergeAlts(self,rows):
        '''merges alt calls for multiple calls at same pos for one sample'''
        r = sorted(rows.ref.values,key=len)[-1]
        a = rows[rows.ref!=r].alt.iloc[0]
        if len(a) == 1: #sub
            alt = a + r[1:]
        else: #ins
            alt = a
        return pd.Series({'ref':r,'alt':alt})

    def new_record(self,loc,data):
        #alts df
        altFunc     = self._getAlt(*loc)
        alts        = self._makeAlts(data,altFunc)
        alleles     = self._getAlleles(alts)
        
        #variant record
        vrec = self.vcfFile.new_record()
        vrec.chrom,vrec.pos = loc
        vrec.id      = data.ID.sort_values(na_position='last').iloc[0] if 'ID' in data.columns else '.' 
        vrec.alleles = alleles
        vrec.qual    = cfg.vcf['defaultQual']
        #TODO make filters more specific
        #eg if all consensus variants are from pbaa-fail clusters, mark variant as pbaafail
        vrec.filter.add('PASS')
        vrec.info['NS'] = len(self.samples) #could be improved to better reflect coverage
        vrec.info['AF'] = list((alts.groupby('alt').size() / len(alts)).reindex(alleles[1:]).values)
        
        #samples
        #for sample,srec in vrec.samples.items():
        for sample in data[self.sampleCol].unique():
            srec        = vrec.samples[sample]
            calls       = data.query(f'{self.sampleCol}==@sample').sort_values(['guide','cluster'])
            merged = False
            if self.mergeVars:
                grouped = calls.groupby('VAR',as_index=False)
                if (grouped.size()['size'] > 1).any(): #some merging
                    calls = grouped.apply(self.mergeVar)\
                                   .reset_index(level=0,drop=True)
                    merged = True
            
            alleleMeta  = self.alleles.loc[calls.index.get_level_values('uuid')]

            suppCount   = calls[cfg.database['supportField']].map(Counter).sum()
            countMap    = {altFunc(v).get('alt',alleles[0]):suppCount[v]
                           for v in ['.'] + calls.reindex(alts.index.get_level_values('uuid'),level=0).VAR.to_list()
                           if v in suppCount}

            #Genotype
            srec['GT']  = list(alts.reset_index(['CHR','POS'],drop=True)\
                                   .reindex(calls[~calls.index.duplicated(keep='first')].index.get_level_values('uuid'))\
                                   .fillna(alleles[0])\
                                   .alt.map(alleles.index))
            srec.phased = not merged 
            #total depth 
            srec['DP'] = int(calls.numreads.astype(int).sum())
            #pbaa avg_qual
            srec['AQ'] = list(alleleMeta.avg_quality) if not merged else None
            #allele depth
            try:
                #srec['AD'] = list(calls.nsupport)
                srec['AD'] = [countMap.get(a,0) for a in alleles]
            except AttributeError:
                print('WARNING: No HiFi support info.  Setting AD to consensus depth')
                srec['AD'] = list(calls.numreads)   
            #variant frequencies
            srec['VAF'] = list(alleleMeta.cluster_freq) if not merged else None
            #target guide
            srec['TG'] = list(calls.guide)
            #cluster number 
            srec['HP'] = list(map(int,alleleMeta.cluster)) if not merged else None
            #pbaa diversity
            srec['DV'] = list(alleleMeta.diversity) if not merged else None
            #pbaa uchime score
            srec['CH'] = list(alleleMeta.uchime_score) if not merged else None
        return vrec

    def mergeVar(self,dat):
        support = cfg.database['supportField']
        f = itemgetter(0)
        d = dat.reset_index()
        firsts = {c:f(d[c]) for c in ['uuid','CHR','POS','VAR']}
        sums   = {c:d[c].astype(int).sum() for c in ['numreads','nsupport']}
        sums[support] = dict(d[support].map(Counter).sum())
        sums['guide'] = ';'.join(d.guide.unique())
        firsts.update(sums)
        return pd.DataFrame([firsts]).set_index(['uuid','CHR','POS'])
         
    def _getAlt(self,ctg,pos):
        ref   = self.reference.fetch(ctg,pos-1,pos).upper()
        names = ['ref','alt']
        def alt(vnt):
            if vnt == '.': #refcall
                return {}
            elif vnt.startswith('*'):
                vals = ref,vnt[-1].upper()
            elif vnt.startswith('-'):
                a=vnt[1:].upper()
                vals = ref+a,ref
            elif vnt.startswith('+'):
                vals = ref,ref+vnt[1:].upper()
            elif vnt.startswith('~'):
                size  = int(vnt[1:-3])
                r_del = self.reference.fetch(ctg,pos-1,pos+size)
                vals  = r_del,r_del[0]
            else:
                raise ValueError(f'unk var {vnt}')
            return dict(zip(names,vals))
        return alt

class Pbaa2Vcf_Error(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='pbaa2vcf.py', description='Generate vcf file from pbAA results')
    parser.add_argument('alleles', metavar='alleles', type=str,
                    help='Alleles csv from consensusVariants.py OR name of alleles table in database.')
    parser.add_argument('variants', metavar='variants', type=str,
                    help='Variants csv from consensusVariants.py OR name of variants table in database')
    parser.add_argument('reference', metavar='reference', type=str,
                    help='Reference fasta (must be indexed)')
    parser.add_argument('-s','--sampleCol', dest='sampleCol', type=str, choices=['barcode','bioSample','runName'], default='barcode',
                    help=f'Column in allelesCsv to use for grouping samples. Default barcode')
    parser.add_argument('-o','--out', dest='out', type=str, default=None,
                    help=f'Output file. Default stdout')
    parser.add_argument('-f','--minFreq', dest='minFreq', type=float, default=cfg.vcf['minFreq'],
                    help=f'Ignore failed clusters below minFreq. Default {DEFAULTMINFREQ}')
    parser.add_argument('-p','--passOnly', dest='passOnly', action='store_true', default=False,
                    help=f'Ignore pbAA failed clusters. Default use minFreq')
    parser.add_argument('-z','--compress', dest='compress', action='store_true', default=False,
                    help=f'Export bcf compressed file [NOT IMPLEMENTED]. Default export uncompressed vcf')
    parser.add_argument('-m','--merge', dest='merge', action='store_true', default=False,
                    help=f'Merge same variants within sample [beta]. Default report all')
    parser.add_argument('-d','--database', dest='database', type=str, default=None,
                    help=f'Sqlite database file.  Default None (use csv inputs)')
    parser.add_argument('-q','--query', dest='query', type=str, default=None,
                    help=f'Custom sql query for input from database. Must return list of uuid values to select on alleles/variant table.  Default None (all records)')

    try:
        vcf = main(parser)
    except Pbaa2Vcf_Error as e:
        print(f'\nERROR: {e}\n')
        sys.exit(1)


