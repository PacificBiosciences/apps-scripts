rule align_consensus:
    input:
        cons=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        ref=config['reference'],
    output:
        bam=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.aligned_maskD7.bam',
        idx=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.aligned_maskD7.bam.bai',
    threads: 
        1
    benchmark:
        f'batches/{batch}/benchmarks/pbmm2/align_consensus.{{sample}}.log'
    conda:
        'envs/pbmm2.yaml'
    shell:
        '''
        pbmm2 align -j {threads} \
                    --preset hifi \
                    --sort \
                    {input.ref} {input.cons} {output.bam}
        '''

rule align_HiFi:
    input:
        bam=f'batches/{batch}/demux/demultiplex.{{barcode}}.bam',
        ref=config['reference'],
    output:
        bam=temp(f'batches/{batch}/aligned/demultiplex.{{barcode}}.aligned_maskd7.bam'),
        idx=temp(f'batches/{batch}/aligned/demultiplex.{{barcode}}.aligned_maskd7.bam.bai'),
    threads: 
        8
    benchmark:
        f'batches/{batch}/benchmarks/pbmm2/align_hifi.{{barcode}}.log'
    conda:
        'envs/pbmm2.yaml'
    shell:
        '''
        pbmm2 align -j {threads} \
                    --preset hifi \
                    --sort \
                    {input.ref} {input.bam} {output.bam}
        '''
