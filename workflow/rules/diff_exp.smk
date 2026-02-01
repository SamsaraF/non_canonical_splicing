# -*- coding=utf-8 -*-

rule build_matrix:
    input:
        expand(
            'results/{{dataset}}/star_align/{sample_name}/ReadsPerGene.out.tab',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/gene_counts_all.tsv'
    params:
        samples = samples['sample_name'].tolist(),
        strandness = lambda wildcards: config[wildcards.dataset]['lib']['strandness']
    script: '../scripts/build_count_matrix.py'


rule deseq2_qc:
    input:
        'results/{dataset}/gene_counts_all.tsv'
    output:
        cluster = 'results/{dataset}/deseq2/cluster.pdf',
        pca = 'results/{dataset}/deseq2/pca.pdf'
    params:
        samples = lambda wildcards: config[wildcards.dataset]['samples'],
        design = lambda wildcards: config[wildcards.dataset]['diff_qc']['design']
    script: '../scripts/deseq2_qc.R'


rule deseq2_diffexp:
    input:
        'results/{dataset}/gene_counts_all.tsv'
    output:
        norm = 'results/{dataset}/deseq2/diffexp_norm.{contrast}.{subset}.tsv',
        ma_plot = 'results/{dataset}/deseq2/diffexp_ma.{contrast}.{subset}.pdf',
        volcano = 'results/{dataset}/deseq2/diffexp_volcano.{contrast}.{subset}.pdf',
    params:
        samples = lambda wildcards: config[wildcards.dataset]['samples'],
        id_to_symbol = lambda wildcards: config[wildcards.dataset]['references']['id_to_symbol'],
        filter_by = lambda wildcards: config[wildcards.dataset]['diff_exp']['filter_by'][wildcards.subset],
        filter_value = '{subset}',
        contrast = lambda wildcards: config[wildcards.dataset]['diff_exp']['contrasts'][wildcards.contrast]
    script: '../scripts/deseq2_diffexp.R'
