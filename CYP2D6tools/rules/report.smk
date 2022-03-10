rule star_typing_all:
    input:
        cons=expand(f'batches/{batch}/{{sample}}/pbaa_{{status}}_cluster_sequences.fasta', \
                     sample=_get_demuxed_samples, status=['passed','failed']),
        biosamples=config['biosamples'],
    output:
        f'batches/{batch}/star_typing_details.csv',
        f'batches/{batch}/star_typing_diplotypes.csv',
    params:
        runname=config['runname'],
        prefix=f'batches/{batch}/star_typing',
    threads: 
        1
    conda:
        'envs/python.yaml'
    shell:
        '''
        python pbCYP2D6typer2.py -H \
                                 -r {params.runname} \
                                 -s {input..biosamples} \
                                 -p {params.prefix} \
                                 {input.cons}
        '''

rule consensus_vcf:
    input:
        cons=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        reads=_get_fastq,
    output:
        f'batches/{batch}/{{sample}}/consensus.vcf',
    params:
        prefix=f'batches/{batch}/{{sample}}/consensus'
    threads: 
        1
    conda:
        'envs/python.yaml'
    shell:
        '''
        python pbCYP2D6typer2.py -H \
                                 --vcf \
                                 -r {params.runname} \
                                 -s {input..biosamples} \
                                 --hifiSupport {input.reads} \
                                 --read_info {input.info} \
                                 -p {params.prefix} \
                                 {input}
        '''
