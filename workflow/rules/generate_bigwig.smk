# -*- coding=utf-8 -*-

rule link_bam:
    input: get_sorted_bam
    output: 'results/{dataset}/bam/{condition_1}.{condition_2}.{replicate}.bam'
    shell: 'ln -s $(realpath {input}) {output}'


rule index_bam:
    input: rules.link_bam.output
    output: 'results/{dataset}/bam/{condition_1}.{condition_2}.{replicate}.bam.bai'
    threads: 4
    shell: 'samtools index -@ {threads} {input}'


rule generate_bigwig:
    input: 
        bam = rules.link_bam.output,
        bai = rules.index_bam.output
    output: 'results/{dataset}/bigwig/{condition_1}.{condition_2}.{replicate}.bw'
    threads: 8
    params:
        bin_size = 1,
        normalize = 'RPKM'
    shell: 'bamCoverage -b {input.bam} -o {output} -bs {params.bin_size} --normalizeUsing {params.normalize} -ignore chrM -p {threads}'
