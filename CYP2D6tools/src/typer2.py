import numpy as np
import pandas as pd
import sqlalchemy as sqla
from sqlalchemy.orm import Session
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

class CypTyper2:
    '''Updated python/pandas star typer.  
       Database is only used for storing metadata (optionally for upload)
    '''
    
    starVariantTables = ['coreSNV','starSNV','starSNVex']

    def __enter__(self):
        return self    

    def __exit__(self,exc_type,exc_value,traceback):
        pass

    def __init__(self,alleles,variants,
                 dbFile,roi='CYP2D6_gene',
                 summarize=True,dbStore=False,
                 cyp2d7Caller=None):
        self.alleles  = alleles
        self.variants = variants
        self.roi_name = roi
        self.dbFile   = dbFile
        self.db       = self._openDB(dbFile)
        self.dbStore  = dbStore
        self.d7caller = cyp2d7Caller
        self.star     = self.getStarVariants()
        self.counts   = self.countVariants()
        self.matches  = self.getMatchSets()
        if summarize:
            self.summary  = self.summarizeResults()
            if self.d7caller is not None:
                self._refineHybrids()
            self.DUP      = self.groupSVhap('DUP').sort_index()
            dupMembers    = [] if self.DUP.empty else self.DUP.members.sum()
            self.HYB      = self.groupSVhap('HYB', exclude=dupMembers).sort_index()
            self.SV       = pd.concat([self.DUP,self.HYB]).sort_index()
            svMembers    = [] if self.SV.empty else self.SV.members.sum()
            self.singles  = self.summary[(~self.summary.uuid.isin(svMembers))
                                         & (self.summary.status == 'passed')].sort_index()
            self.diplotypes = self.makeDiplotypes()
        if dbStore:
            #TODO: allow upload to different db
            self._uploadDB(dbFile)
        
    @property
    def roi(self):
        '''Return (start,stop) POS for roi'''
        return self.star['roi'].query('name==@self.roi_name')[['start','stop']].iloc[0]
    
    def _openDB(self,dbFile):
        '''Open db connection'''
        return sqla.create_engine(f'sqlite:///{dbFile}', echo=False)
    
    def _uploadDB(self,destDBfile):
        '''Import alleles and variants to typing database'''
        with sqla.create_engine(f'sqlite:///{destDBfile}', echo=False) as engine:
            tables   = [cfg.database[t] for t in ['alleleTable','variantTable']]
            support  = cfg.database['supportField'] 
            maxTries = cfg.database['maxTries']
            for pydf,sqldf in zip([self.alleles,self.variants],tables):
                pbaa = pydf.reset_index()\
                           .reindex(columns=list(cfg.tableMap[sqldf].values()))
                pbaa.columns = list(cfg.tableMap[sqldf].keys())
                #make sure support dtyep is "string"
                if support in pbaa.columns:
                    pbaa[support] = pbaa[support].astype(str)
                tries = 0
                while tries < maxTries:
                    try:
                        pbaa.set_index('uuid').to_sql(sqldf, con=self.db, if_exists='append')
                        break
                    #except sqlite3.OperationalError as e:
                    except sqla.exc.OperationalError as e:
                        tries += 1
                        print(f'WARNING: sqlite3 import error try #{tries} of {maxTries}')
                        if tries == maxTries:
                            raise Typer_Error(f'Unable to import {pydf.source.unique()} to {self.dbFile}')
                        else:
                            time.sleep(1)
            return None
            
    
    def getStarVariants(self):
        '''Gather star-allele data from the db'''
        print(f'Reading from {self.dbFile}')
        tables = {'coreSNV'    : 'STAR_core_snv_roi',
                  'coreIMPACT' : 'CYP2D6_impacts',
                  'starVar'    : 'STAR_variants_roi',
                  'cyp2d6Exons': 'CYP2D6_exons',
                  'cyp2d7Exons': 'CYP2D7_exons',
                  'starAlleles': 'STAR_alleles',
                  'roi'        : 'ROI'}
        #dict of name : dataframe
        star = { name : pd.read_sql_table(table,con=self.db)
                 for name,table in tables.items() }
        #remove rows from core with no variant
        star['coreSNV'] = star['coreSNV'].query('VAR.notnull()',engine='python')
        #merge allele meta with snv for all gDNA variants
        starSNV = pd.merge(star['starVar'],star['starAlleles'],on='uuid')\
                    .query('VAR.notnull()',engine='python')\
                    .drop('uuid',axis=1)
        star['starSNV'] = starSNV
        #mask intronic variants
        isExonic = lambda pos: (len(star['cyp2d6Exons'].query('start <= @pos <= stop')) > 0)
        star['starSNVex'] = starSNV[starSNV.POS.map(isExonic)]
        return star
    
    
    def _getCol(self,name):
        '''utility func to map table type to unique column'''
        return 'CoreAllele' if 'core' in name else 'alleleName'
    
    def countVariants(self):
        '''counts of star variants over ROI.  Use this to filter amplicons without all expected calls'''
        print('Counting variants')
        #star variants
        counts = { name : table.groupby(self._getCol(name)).size().rename('refVarCount')
                   for name,table in self.star.items()
                   if name in self.starVariantTables }
        #query samples
        start,stop = self.roi
        counts['samples'] = self.variants.query('@start <= POS <= @stop')\
                                         .groupby('uuid').size()\
                                         .rename('sampleVarCount')
        return counts
        
    def getMatchSets(self):
        '''call all matching star alleles core/gdna/cdna'''
        print('Generating match sets')
        return { name : self.matchVarSets(name,ref,self.variants.reset_index())
                 for name,ref in self.star.items()
                 if name in self.starVariantTables}
    
    def matchVarSets(self,refname,reference,query):
        '''filter calls without the right number of vars and sort by score/nvar'''
        refColumn = self._getCol(refname)
        grpCols   = ['uuid',refColumn]
        isCore    = 'core' in refname
        #cross join of ref and query
        sharedVars   = pd.merge(reference,query,
                                how='inner',on=['CHR','POS','VAR'])
        #compare variant set count, keep only candidates with all ref variants over roi
        sharedCounts = sharedVars.groupby(grpCols).size()\
                                 .rename('varCountShared')\
                                 .reset_index('uuid')
        refCounts    = self.counts[refname]
        matchedCount = pd.merge(sharedCounts,refCounts,
                                on=refColumn)\
                         .set_index('uuid',append=True)\
                         .query('varCountShared == refVarCount')
        if isCore: 
            matchedCount.reset_index(refColumn,inplace=True)
        else:#need to match counts on sample as well for non-core matches
            matchedCount = matchedCount.reset_index(refColumn)\
                                       .join(self.counts['samples'])\
                                       .query('varCountShared == sampleVarCount')
        #get core/major allele impacts
        cypMajor = matchedCount[refColumn].str.split('.').str[0].values
        impacts  = self.star['coreIMPACT']\
                       .set_index('alleleMajor')\
                       .reindex(cypMajor)\
                       .reset_index()
        impacts.index = matchedCount.index
        
        matches = pd.concat([matchedCount,impacts],axis=1)\
                    .eval('score = impactVal + priorityMod')\
                    .groupby('uuid').apply(self.sortMatches(refColumn))
        
        return pd.Series(dtype='object') if matches.empty else matches
        
    def parseCA(self,ca):
        '''Parse names of the form "CYP2D6_3" --> "*3" '''
        return f'*{ca.split("_")[-1]}'

    def sortMatches(self,colname):
        '''returns the match sorting function'''
        def sFunc(matches):
            'select highest score followed by most shared vars'
            sortMatches = matches.sort_values(['score','varCountShared'],ascending=False)
            return ';'.join(map(self.parseCA,sortMatches[colname]))
        return sFunc    
        
    def summarizeResults(self):
        '''merge calls, select best core, label SV, label *1'''
        print('Summarizing results')
        allMatches = pd.concat({s:self.matches[s] for s in self.starVariantTables},axis=1,names=['uuid'])
        withMeta   = pd.concat([allMatches,self.alleles],axis=1)
        withMeta.index.name = 'uuid'
        withMeta['SV'] = None
        #label *5
        withDels   = self.label5Deletion(withMeta)
        withSV     = self.labelSV(withDels)
        finalCalls = self.label1(withSV)
        finalCalls['Type'] = finalCalls.coreSNV.str.split(';').str[0]
        
        #select output cols and fix esthetics
        colNames = OrderedDict(bioSample='bioSample',
                               barcode='barcode',
                               uuid='uuid',
                               cluster='cluster',
                               numreads='numreads',
                               length='seqLength',
                               hasSpacer='hasSpacer',
                               Type='Type',
                               SV='SV',
                               coreSNV='core_all',
                               starSNV='gDNA_all',
                               starSNVex='cDNA_all',
                               clusterStatus='status')
        return finalCalls.reset_index()\
                         .rename(columns=colNames)[list(colNames.values())]\
                         .set_index(['bioSample','barcode'])\
                         .sort_index()
    
    def label5Deletion(self,table):
        '''find deletion alleles in set and label.
           TODO: add test of gapped mapping of amplicon
        '''
        delstar5 = table.query('coreSNV.isnull() \
                                and 4000 < length < 6000 \
                                and alnStart > 42123000 \
                                and alnStop < 42132000',engine='python').index
        table.loc[delstar5,'coreSNV'] = '*5'
        table.loc[delstar5,'SV'] = 'DEL'
        return table
    
    def labelSV(self,table):
        '''heuristics to identify dups and hybrids'''
        dups = table.query('hasSpacer == 1 \
                            and 8200 < length < 8800',engine='python').index
        hybs = table.query('hasSpacer == 1 \
                            and length > 9000',engine='python').index
        table.loc[dups,'SV'] = 'DUP'
        table.loc[hybs,'SV'] = 'HYB'
        return table
    
    def label1(self,table):
        '''remaining alleles with no core call that span the gene are labeled *1'''
        refcalls = table.query('coreSNV.isnull() \
                                and 7800 < length < 9000 \
                                and alnStart < 42126400 \
                                and alnStop > 42130850',engine='python').index
        table.loc[refcalls,'coreSNV'] = '*1'
        return table
    
    def getJaccardSimilarity(self,uuid1,uuid2):
        '''jaccard sim of variants for two alleles'''
        try:
            union = pd.concat([self.variants.loc[uuid1],
                               self.variants.loc[uuid2]],
                               axis=1)
        except KeyError as e: #some uuid has no variants
            return 0
        if union.empty:
            return 0
        samePos   = union.dropna()
        intersect = samePos[(samePos.iloc[:,0] == samePos.iloc[:,1])].iloc[:,0]
        return len(intersect) / len(union)
    
    def groupSVhap(self,kind='DUP',exclude=None):
        '''Generate haplotype call for multi-allele haplotypes, eg "*10 + *36" or "*2 x2"'''
        if kind == 'DUP':
            sep = ' x'
            svf = lambda d: len(d) + 1 
            on  = ['bioSample','barcode','Type']
        elif kind == 'HYB':
            sep = ' + '
            svf = lambda d: sep.join(d.Type_HYB)
            on  = ['bioSample','barcode']
        else:
            raise ValueError
            
        def makename(grp):
            '''generate name, given type and uuids'''
            members = tuple(set(list(grp.uuid.values) + list(grp[f'uuid_{kind}'].values)))
            allele  = pd.Series({'members'  : members,
                                 'haplotype': f'{grp.Type.iloc[0]}{sep}{svf(grp)}'})
            return allele 
    
        svMatch = pd.merge(self.summary.query('SV.isnull() and status=="passed"',engine='python'),
                           self.summary.query('SV==@kind and status=="passed"'),
                           on=on,
                           suffixes=[None,f'_{kind}'])
        if exclude is not None:
            svMatch = svMatch[~svMatch.uuid.isin(exclude)]

        if svMatch.empty:
            for col in ['jaccardSim','members','haplotype']:
                svMatch[col] = None
            return svMatch
        else:
            svMatch['jaccardSim'] = svMatch.apply(lambda row: self.getJaccardSimilarity(row.uuid,row[f'uuid_{kind}']),axis=1)
            svMatch.sort_values(['bioSample','barcode',f'uuid_{kind}','jaccardSim','numreads'],
                                ascending=False,inplace=True)
        
            return svMatch[~svMatch.duplicated(f'uuid_{kind}',keep='first')]\
                          .groupby('uuid').apply(makename)\
                          .join(self.alleles[on[:2]])\
                          .set_index(on[:2]).sort_index()
    
    def makeDiplotypes(self):
        '''Generate diplotypes'''
        def munger(idx):
            single = self.singles.Type.get(idx,'')
            if type(single) == pd.Series:
                # slightly hacky "dropna" to remove odd fragments with no type 
                # these will still be in the detailed summary
                single = list(single.dropna().values)
            else:
                single = [single]
            sv = self.SV.haplotype.get(idx,'')
            if type(sv) == pd.Series:
                sv = list(sv.values)
            else:
                sv = [sv]
            haplotypes = pd.Series(single+sv)
            haplotypes = haplotypes[haplotypes != '']
            
            sortedHap  = pd.concat([haplotypes,
                                    haplotypes.str.split(' ')\
                                              .str[0].str[1:]\
                                              .astype(int).rename('coreNumber')],
                                   axis=1)\
                           .sort_values('coreNumber')\
                           .iloc[:,0].values
            
            return ' / '.join(sortedHap)

        index = self.singles.index.union(self.SV.index).unique()
        return pd.Series(map(munger,index),
                         index=index,
                         name='diplotype')

    def _refineHybrids(self):
        #replicate variant calling while masking cyp2d6
        def getHyb(uuid):
            try:
                varnt6   = self.variants.loc[uuid].reset_index('POS')
                varnt7   = self.d7caller.variants.loc[uuid].reset_index('POS')
                hycaller = hybridCaller(varnt6,varnt7,
                                        self.star['cyp2d6Exons'],
                                        self.star['cyp2d7Exons'],
                                        self.alleles.loc[uuid],
                                        names=['CYP2D6','CYP2D7'],
                                        countBases=True,
                                        nwindows=35)
                return hycaller
            except KeyError:
                return None  

        self.summary['hybrid_exons'] = None
        
        for uuid in self.summary.query('SV=="HYB"').uuid:
            hycaller = getHyb(uuid)
            if hycaller:
                mask = self.summary.uuid == uuid
                self.summary.loc[mask, 'Type'] = hycaller.starAllele
                self.summary.loc[mask, 'hybrid_exons'] = f'{hycaller.exonSummary.to_dict()}'
                #TODO would be nice to print hybrid plot here

class hybridCaller:
    color = {'A':'r','B':'b'}

    def __init__(self,variants_A,variants_B,
                      exonTable_A,exonTable_B,
                      alleleMeta,
                      names=None,
                      sameStrand=True,
                      countBases=True,
                      nwindows=500
                 ):
        '''*_A        : alignments and start/stop for primary gene
           *_B        : alignments and start/stop for secondary gene
           sameStrand : same strand (True) or reverse the secondary (False)
        '''
        self.keys     = list('AB')
        self.var      = {'A': variants_A, 'B': variants_B}
        self.exons    = {'A': exonTable_A, 'B': exonTable_B}
        self.names    = {'A':names[0] if names else 'A', 'B':names[1] if names else 'B'}
        self.meta     = alleleMeta
        #self.same     = sameStrand
        self.countFunc = self._countVarBases if countBases else self._countVarPos
        self.win      = self._getWindows(nwindows)
        self.predPath,self.switch = self._get_predPath()
        self.exonSim  = self._get_simExons()
        
    @property    
    def roi(self):
        return pd.DataFrame({ k : self.exons[k].apply({'start':min,'stop':max})
                              for k in self.keys })
    
    def _getWindows(self,nwindows):
        return { k : pd.interval_range(self.roi.loc['start',k],self.roi.loc['stop',k],nwindows)
                 for k in self.keys }
    
    @property
    def diffs(self):
        intervals = [self.win[t] for t in self.keys]
        return pd.concat({k: pd.Series({ intv.get_loc(i) : n 
                                         for i,n in self.countFunc(self.var[k],intv).items() 
                                       })
                           for k,intv in zip(self.keys,intervals)},
                          axis=1)
    
    def _baseDiffs(self,vnt):
        if vnt.startswith('*'):
            return 1
        else:
            return len(vnt) - 1
    
    def _countVarPos(self,variants,interval):
        return pd.cut(variants.POS,interval).dropna().value_counts()

    def _countVarBases(self,variants,interval):
        baseCounts = variants.VAR.map(self._baseDiffs)
        return baseCounts.groupby(pd.cut(variants.POS,interval)).sum()
    
    @property
    def fwdDiff(self):
        return self.diffs.cumsum().rename(columns=self.names)
    
    @property
    def revDiff(self):
        return self.diffs[::-1].cumsum()[::-1].rename(columns=self.names)
    
    @property
    def pathDiff(self):
        names = list(self.names.values())
        return pd.concat({(fn,rn):pd.concat([self.fwdDiff[fn],self.revDiff[rn]],axis=1).min(axis=1)
                          for fn in names for rn in names},axis=1)
    
    def _get_predPath(self):
        '''new selection algo'''
        path   = pd.concat([self.fwdDiff,self.revDiff],axis=1).idxmin(axis=1)
        switch = path[(path != path[0]).cumsum() == 1]
        return path,switch
        #old version
        #idx    = self.pathDiff.sum().idxmin()
        #joined = pd.concat([self.fwdDiff[idx[0]],self.revDiff[idx[1]]],axis=1)
        #path   = joined.idxmin(axis=1)
        #switch = path[(path != path[0]).cumsum() == 1]
        #return path,switch
    
    @property
    def mostSim(self):
        'Filling in ambiguous regions with predicted path, otherwise whichever has fewest vars'
        return self.predPath.where(self.diffs.diff().B == 0,
                                   self.diffs.idxmin(axis=1).map(self.names))
    
    def _get_simExons(self,tol=2):
        res = {}
        for k in self.keys:
            varPos     = self.var[k].POS.values
            exonStarts = (self.exons[k].start - tol).values[:,None]
            exonStops  = (self.exons[k].stop + tol).values[:,None]
            counts     = ((varPos >= exonStarts) & (varPos <= exonStops)).sum(axis=1)
            index      = self.exons[k].exon.astype(int)
            res[self.names[k]] = pd.Series(counts,index)
        df = pd.DataFrame(res)
        out = df.idxmin(axis=1)
        out[(df.diff(axis=1).iloc[:,-1] == 0)] = 'ambiguous' #exons same-ish
        return out
    
    @property
    def starAllele(self):
        'Temporary Heuristc'
        loci = self.predPath.unique()
        if len(loci) == 1:
            if self.meta.hasSpacer or self.meta.length > 8500:
                return f'Novel {loci[0]} hybrid?'
            else:
                return f'{loci[0]} non-Hybrid'
        if self.exonSummary['CYP2D7'] == (9,): #just exon9 is 2D7
            return '*36'
        if 9 in self.exonSummary['CYP2D6']: # 2D7-2D6 (*13)
            return '*13'
        if 1 in self.exonSummary['CYP2D6']: # *68 or *83
            if 8 in self.exonSummary['CYP2D6']:
                return '*68'
            else:
                return '*83'
        else: #unknown hyb
            return 'unknown-hybrid'
    
    def _get_breakRange(self,tolerance=3):
        #for non-hybrids, return empty list
        if self.switch.empty:
            return []
        else:
            breakRange  = [self.switch.index[0]]
        #guessing what these two values are
        #fromLocus = 1
        #toLocus = 2
        #if fromLocus != toLocus:
            #check backwards
            #backChange  = joined.diff().fillna(joined).abs().min(axis=1)
            backChange  = self.revDiff.diff(axis=1).iloc[:,1].diff().abs()
            region      = self.switch.index[0] - 1
            while backChange[region] <= tolerance:
                breakRange.append(region)
                region -= 1
                if region < 0: break
            #check forwards
            fwdChange = self.fwdDiff.diff(axis=1).iloc[:,1].diff().abs()
            region = self.switch.index[0] + 1
            while fwdChange[region] <= tolerance:
                breakRange.append(region)
                region += 1
                if region == len(fwdChange): break
        return sorted(breakRange)
    
    @property
    def exonSummary(self):
        return self.exonSim.groupby(self.exonSim).apply(lambda d: tuple(d.index.sort_values()))
    
    def _invertIndex(self,s):
        return s[::-1].reset_index(drop=True)
    
    def plotHybrid(self,size=(14,8),scale=0.5):
        shift = 10 #esthetics only
        color = [self.color[k] for k in self.keys]
        f,ax = plt.subplots(figsize=size)
        #plot forward/reverse accumulation of vars
        cols = self.fwdDiff.columns
        fwdCol,revCol = [dict(zip(cols,map(f'{{}}_{d}'.format,cols))) for d in ['fwd','rev']]
        ((self.fwdDiff.rename(columns=fwdCol) + shift)*scale).plot(style='->',ax=ax,color=color,markersize=5)
        ((self.revDiff.rename(columns=revCol) + shift)*scale).plot(style='<:',ax=ax,color=color,markersize=5)
        
        breakRange = self._get_breakRange()
        
        #exons
        intv   = self.win['A']
        cMap   = {self.names[k] : self.color[k] for k in self.keys}
        cMap['ambiguous'] = 'm'
        colors = self.exonSim.map(cMap)
        ex_mid = pd.cut(self.exons['A'].eval('start + (stop-start)/2'),intv).map(intv.get_loc)
        ex_num = self.exons['A'].exon
        ax.vlines(ex_mid,-6,-12,colors=colors)
        for exnum,expos,color in zip(ex_num,ex_mid,colors):
            ax.annotate(f'e{exnum}',(expos,-15),color=color)
            
        #shade regions of rapid differentiation
        deltaRegions = (self.fwdDiff.diff(axis=1).diff().shift(-1).abs().iloc[:,-1] > 15) * (self.fwdDiff.max().max() + shift)
        deltaRegions.plot(kind='bar',width=1,color='k',alpha=0.05,ax=ax)
        
        #most similar over window
        kwargs = dict(kind='bar',width=1,alpha=1,bottom=-2)
        for k,gene in self.names.items():
            ((self.mostSim == gene).astype(int) * -3).plot(color=self.color[k],**kwargs)
        #add purple for ambiguous similarity
        ((self.diffs.diff().B == 0) * -3).plot(color='m',**kwargs)
        ax.annotate('Closest Match',(self.mostSim.index.max()+.5,-8.1))
        #predicted hybrid
        kwargs.pop('bottom')
        for k,gene in self.names.items():
            ((self.predPath == gene).astype(int) * 3).plot(color=self.color[k],**kwargs)
        #add ambiguous break
        (self.predPath.reset_index()['index'].isin(breakRange) * 3).plot(color='m',**kwargs)
        ax.annotate('Hybrid Prediction',(self.predPath.index.max()+.5,2.1))
        
        #add patches to legend
        handles, labels = ax.get_legend_handles_labels()
        d6    = mpatches.Patch(color='r', label='2D6-like')
        d7    = mpatches.Patch(color='b', label='2D7-like')
        ambig = mpatches.Patch(color='m', label='Ambiguous')
        delta = mpatches.Patch(color='k', ec=None, alpha=0.05, label='Rapid Differentiation')
        newhandles = handles[:4] + [delta,d6,d7,ambig] 
        plt.legend(handles=newhandles)#, loc='upper center')
        
        #axes labels
        ax.set_ylabel('Cumulative variants')
        ax.set_xlabel('Gene Region')
        #ax.annotate('Hybrid Prediction',(hybPrediction.index.max()+.5,2.1))
        #ax.annotate('Closest Match',(hybPrediction.index.max()+.5,-5.1))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        #remove tick marks
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        
        ax.set_title(f'{self.starAllele}')
