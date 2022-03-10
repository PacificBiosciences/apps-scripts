rule star_typing_all:
    input:
        cons=_agg_consensus,
        biosamples=config['biosamples'],
    output:
        f'batches/{batch}/{config["prefix"]}_detailed_summary.csv',
        f'batches/{batch}/{config["prefix"]}_diplotype_summary.csv',
    params:
        script=config['starscript'],
        runname=batch,
        prefix=f'batches/{batch}/{config["prefix"]}',
        extra='-H',
    threads: 
        1
    conda:
        'envs/python.yaml'
    shell:
        '''
        python {params.script} -r {params.runname} \
                               -s {input.biosamples} \
                               -p {params.prefix} \
                               {params.extra} \
                               {input.cons}
        '''

rule consensus_vcf:
    input:
        cons=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        biosamples=config['biosamples'],
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        reads=_get_fastq,
    output:
        f'batches/{batch}/{{sample}}/consensus_passed.vcf',
    params:
        script=config['starscript'],
        runname=batch,
        prefix=f'batches/{batch}/{{sample}}/consensus_passed',
        extra='-H --vcf --ignoreMissing'
    threads: 
        1
    conda:
        'envs/python.yaml'
    shell:
        '''
        python {params.script} -r {params.runname} \
                               -s {input.biosamples} \
                               --hifiSupport {input.reads} \
                               --read_info {input.info} \
                               -p {params.prefix} \
                               {params.extra} \
                               {input.cons}
        '''
