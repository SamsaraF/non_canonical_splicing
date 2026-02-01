# -*- coding=utf-8 -*-

def infer_adapter(wildcards):
    adapter = config[wildcards.dataset]['lib'].get('adapter', 'none')
    params_map = {
        'illumina': '--illumina',
        'nextera': '--nextera',
        'small_rna': '--small_rna'
    }
    if adapter.lower() == 'auto':
        return ''
    elif adapter.lower() == 'none':
        raise ValueError('Adapter trimming is set to "none", which means "trim_adapter" rule should be skipped, but was called unexpectedly.')
    elif adapter.lower() in params_map:
        return params_map[adapter.lower()]
    else:
        raise ValueError(f'Unknown adapter type: {adapter}. Supported types are: {list(params_map.keys())}.')


rule trim_adapter:
    input: unpack(get_raw_fastq)
    output:
        r1 = 'results/{dataset}/cleaned_fastq/{sample_name}_val_1.fq.gz',
        r2 = 'results/{dataset}/cleaned_fastq/{sample_name}_val_2.fq.gz'
    params:
        adapter = lambda wildcards: infer_adapter(wildcards),
        out_dir = subpath(output.r1, parent=True)
    threads: 4
    shell:
        'trim_galore -j {threads} --paired {params.adapter} --stringency 3 -q 25 '
        '-o {params.out_dir} --basename {wildcards.sample_name} {input.r1} {input.r2}'
    