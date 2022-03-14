checkpoint demux:
    input:
        ccs=config['hifireads'],
        barcodes=config['barcodeFa'],
        biosamples=config['biosamples']
    output:
        directory(f'batches/{batch}/demux')
    threads:
        24
    benchmark:
        f'batches/{batch}/benchmarks/demux/lima.log'
    conda:
        'envs/lima.yaml'
    shell:
        '''
        mkdir -p {output}
        lima --hifi-preset ASYMMETRIC \
             --split-named \
             -j {threads} \
             --log-level INFO \
             --biosample-csv {input.biosamples} \
             {input.ccs} {input.barcodes} {output}/demultiplex.bam
        '''
