rule prep_HiFi:
    input:
        f'batches/{batch}/demux/demultiplex.{{barcode}}.bam'
    output:
        fq=temp(f'batches/{batch}/fastq/demultiplex.{{barcode}}.fastq'),
        idx=temp(f'batches/{batch}/fastq/demultiplex.{{barcode}}.fastq.fai'),
    threads: 
        1
    benchmark:
        f'batches/{batch}/benchmarks/samtools/bam2fastq.{{barcode}}.log'
    conda:
        'envs/samtools.yaml'
    shell:
        '''
        samtools fastq {input} > {output.fq}
        samtools fqidx {output.fq}
        '''

rule pbaa_cluster:
    input:
        fq=_get_fastq,
        idx=_get_fq_idx,
    output:
        cons1=f'batches/{batch}/{{sample}}/pbaa_passed_cluster_sequences.fasta',
        cons2=f'batches/{batch}/{{sample}}/pbaa_failed_cluster_sequences.fasta',
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        log=f'batches/{batch}/{{sample}}/pbaa.log'
    params:
        prefix=f'batches/{batch}/{{sample}}/pbaa',
        guide=config['guide'],
        loglevel='INFO',
    threads:
        16
    benchmark:
        f'batches/{batch}/benchmarks/pbaa/cluster.{{sample}}.log'
    conda:
        'envs/pbaa.yaml'
    shell:
        '''
        pbaa cluster --min-cluster-frequency 0.05 \
                     --max-uchime-score 10 \
                     --max-reads-per-guide 1000 \
                     --min-var-frequency 0.3 \
                     -j {threads} \
                     --log-file {output.log} \
                     --log-level {params.loglevel} \
                     {params.guide} {input.fq} {params.prefix}
        '''

rule extract_clustered_reads:
    input:
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        bam=_get_align_bam_from_sample,
    output:
        incl=temp(f'batches/{batch}/{{sample}}/clustered_holes.txt'),
        tmp=temp(f'batches/{batch}/{{sample}}/hifi.bam'),
    threads:
        1
    conda:
        'envs/samtools.yaml'
    benchmark:
        f'batches/{batch}/benchmarks/samtoole/extract.{{sample}}.log'
    shell:
        '''
        cut -d' ' -f1 {input.info} > {output.incl}
        samtools view -bhN {output.incl} {input.bam} -o {output.tmp}
        '''
        
rule paint_bam:
    input:
        info=f'batches/{batch}/{{sample}}/pbaa_read_info.txt',
        bam=_get_align_bam_from_sample
    output:
        f'batches/{batch}/{{sample}}/hifi.painted.bam'
    threads: 
        1
    conda:
        'envs/pbaa.yaml'
    benchmark:
        f'batches/{batch}/benchmarks/pbaa/bampaint.{{sample}}.log'
    shell:
        '''
        pbaa bampaint {input.info} {input.bam} {output}
        '''

rule index_bam:
    input:
        f'batches/{batch}/{{sample}}/hifi.painted.bam'
    output:
        f'batches/{batch}/{{sample}}/hifi.painted.bam.bai'
    threads: 
        1
    benchmark:
        f'batches/{batch}/benchmarks/samtools/index_painted.{{sample}}.log'
    conda:
        'envs/samtools.yaml'
    shell: 
        '''
        samtools index {input}
        '''
