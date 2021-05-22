#!/usr/bin/env python3
__version__ = '0.1.0'
import os,sys
from operator import xor
from datetime import datetime
from src.variants import Caller,Aligner
from src.typer import CypTyper
from src import config as cfg

def main(parser):
    args = parser.parse_args()

    if args.vcf and (args.hifiSupport is None or args.read_info is None):
        raise pbCYP2D6typer_Error('VCF export requires both hifiSupport and read_info')

    if xor(args.hifiSupport is None,args.read_info is None):
        raise pbCYP2D6typer_Error('Use both hifiSupport and read_info together')

    reference = getPath('reference')
    spacer    = getPath('spacer') #unique spacer assoc w/D7

    vCaller = Caller(args.consensusFasta,
                     args.runName,
                     reference,
                     spacer,
                     args.preset,
                     args.sampleMap,
                     args.minFrac,
                     args.ignoreMissing,
                     args.read_info,
                     args.hifiSupport,
                     getNow())

    if args.vcf:
        from src.vcf import VcfCreator
        vcf = VcfCreator(f'{args.prefix}.vcf',
                         vCaller.alleles,
                         vCaller.variants,
                         reference,
                         sampleCol='barcode' if args.sampleMap is None else 'bioSample',
                         passOnly=True,
                         dataframe=True) 
        vcf.run()

    if args.writeCSV:
        vCaller.alleles.to_csv(f'{args.prefix}_alleles.csv')
        vCaller.variants.to_csv(f'{args.prefix}_variants.csv')        

    with CypTyper(vCaller.alleles,vCaller.variants,\
                  getPath('database'),args.dbStore) as starTyper:

        for summary,descr in zip(['summaryTable','shortSummary'],['detailed','short']):    
            starTyper.getSummary(cfg.database[summary])\
                     .to_csv(f'{args.prefix}_{descr}_summary.csv',index=False)
        
        #TODO: export all potential matches
        #if args.allCandidateTable:
            #get long report

    return vCaller

def getNow():
    now = datetime.now()
    return now.strftime(cfg.caller['dateFormat'])

def getPath(item):
    scriptPath = os.path.abspath(os.path.dirname(sys.argv[0]))
    return f'{scriptPath}/{cfg.dataPaths[item]}'

class pbCYP2D6typer_Error(Exception):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog='pbCYP2D6typer.py', description='Call CYP2D6 star alleles from pbaa consensus')
    parser.add_argument('consensusFasta', metavar='consensusFasta', nargs='*', type=str,
                    help='Fasta file(s) of pbAA consensus outputs')
    inputp = parser.add_argument_group('Input Options')
    inputp.add_argument('-r','--runName', dest='runName', type=str, default=None, required=True,
                    help=f'Sequencing run ID. Required')
    inputp.add_argument('-s','--sampleMap', dest='sampleMap', type=str, default=None,
                    help=f'CSV file mapping biosample name to barcode. Must have fields {cfg.caller["sampleMapCols"]}. Default None')
    inputp.add_argument('-i','--ignoreMissing', dest='ignoreMissing', action='store_true', default=False,
                    help=f'Do not report empty rows for missing samples from sample map. Default False')
    inputp.add_argument('-f','--minFrac', dest='minFrac', type=float, default=cfg.caller['minFrac'],
                    help=f'Ignore failed clusters below minFrac. Default {cfg.caller["minFrac"]}')
    inputp.add_argument('-P','--preset', dest='preset', choices=list(Aligner.presets.keys()), default=cfg.caller["preset"],
                    help=f'Alignment preset for mappy aligner. Default {cfg.caller["preset"]}')
    inputp.add_argument('--hifiSupport', dest='hifiSupport', type=str, default=None,
                    help=f'Add per-variant read support depth. Path to hifi fastq used for clustering. Requires read_info option to be set. Default None')
    inputp.add_argument('--read_info', dest='read_info', type=str, default=None,
                    help=f'Path to pbAA read_info table for the run. Required for hifi-support option. Default None')
    outputp = parser.add_argument_group('Output Options')
    outputp.add_argument('-p','--prefix', dest='prefix', type=str, default='./cyp2d6typer', required=False,
                    help=f'Output prefix. Default ./cyp2d6typer')
    outputp.add_argument('--writeCSV', dest='writeCSV', action='store_true', default=False,
                    help=f'Write intermediate CSV files. Default False')
    outputp.add_argument('--vcf', dest='vcf', action='store_true', default=False,
                    help=f'Output VCF file. HiFi support and read_info required. Default False')
    outputp.add_argument('--dbStore', dest='dbStore', action='store_true', default=False,
                    help=f'Maintain sample records in local db. Default False')

    try:
        vCaller = main(parser)
    except pbCYP2D6typer_Error as e:
        print(f'\nERROR: {e}\n')
        sys.exit(1)



