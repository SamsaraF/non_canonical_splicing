# -*- coding=utf-8 -*-

rule fastq_dump:
    input: 'data/{dataset}/{accession}.sra'
    output:
        r1 = 'results/{dataset}/fastq/{accession}_1.fastq',
        r2 = 'results/{dataset}/fastq/{accession}_2.fastq'
    threads: 8
    params:
        out_dir = subpath(output.r1, parent=True)
    shell:
        'fasterq-dump -e {threads} --split-3 -O {params.out_dir} {input}'


rule gzip_fastq:
    input: 'results/{dataset}/fastq/{accession}_{read}.fastq'
    output: 'results/{dataset}/fastq/{accession}_{read}.fastq.gz'
    threads: 8
    shell: 'pigz -p {threads} {input}'
