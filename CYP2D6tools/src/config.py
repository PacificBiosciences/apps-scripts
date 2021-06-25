__version__ = '0.1.3'

dataPaths = {'reference': 'src/db/human_GRCh38_no_alt_analysis_set_chr22only_cyp2d7-masked.fasta',
             'database' : 'src/db/CYP2D6.db',
             'spacer'   : 'src/db/spacer.fasta'}

database = {'alleleTable' : 'pbAA_consensus',
            'variantTable': 'SampleVariants',
            'variantAnnot': 'SampleVariants_annotated',
            'summaryTable': 'type_summary',
            'shortSummary': 'type_summary_short',
            'supportField': 'supportReads',
            'maxTries'    : 10}

tableMap = {database['alleleTable']:      {'uuid'               :'uuid',
                                           'runName'            :'runName',
                                           'bioSample'          :'bioSample',
                                           'barcode'            :'barcode',
                                           'guide'              :'guide',
                                           'cluster'            :'cluster',
                                           'numreads'           :'numreads',
                                           'uchime_score'       :'uchime_score',
                                           'uchime_left_parent' :'uchime_left_parent',
                                           'uchime_right_parent':'uchime_right_parent',
                                           'cluster_freq'       :'cluster_freq',
                                           'diversity'          :'diversity',
                                           'avg_quality'        :'avg_quality',
                                           'filters'            :'filters',
                                           'source'             :'source',
                                           'clusterStatus'      :'clusterStatus',
                                           'chrom'              :'chrom',
                                           'alnStart'           :'alnStart',
                                           'alnStop'            :'alnStop',
                                           'hasSpacer'          :'hasSpacer',
                                           'csString'           :'csString',
                                           'length'             :'length',
                                           'status'             :'clusterStatus',
                                           'datetime'           :'datetime'},
            database['variantTable']:     {'uuid'                    :'uuid',
                                           'CHR'                     :'CHR',
                                           'POS'                     :'POS',
                                           'VAR'                     :'VAR',
                                           database['supportField']  :database['supportField']}}

caller = {'namePattern'  : 'sample-(?P<barcode>.*)_guide-(?P<guide>.*)_cluster-(?P<cluster>[0-9]+)_ReadCount-(?P<numreads>[0-9]+)',
          'minFrac'      : 0.01,
          'minLength'    : 2000,
          'preset'       : 'gapstrict',
          'dateFormat'   : '%Y-%m-%d %H:%M:%S',
          'readInfo'     : ['readName',
                            'guide',
                            'orient',
                            'secondBestGuide',
                            'Score',
                            'ScoreParts',
                            'Sample',
                            'length',
                            'avgQuality',
                            'ClusterId',
                            'ClusterSize'],
          'sampleMapCols': ['Barcode','Bio Sample Name']}

typer = {'sortColumns' : ['bioSample','barcode','type','cluster']}

vcf = {'defaultQual' : 200,
       'minFreq'     : 0.05}
