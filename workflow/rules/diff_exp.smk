# -*- coding=utf-8 -*-

rule build_matrix:
    input:
        expand(
            'results/{{dataset}}/star_align<star_method>/{sample_name}/ReadsPerGene.out.tab',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/star_align<star_method>/star_read_counts.all.tsv'
    params:
        samples = samples['sample_name'].tolist(),
        strandness = lambda wildcards: config[wildcards.dataset]['lib']['strandness']
    script: '../scripts/build_count_matrix.py'


FEATURECOUNTS_STRANDNESS_MAP = {
    'un': '0',
    'fr': '1',
    'rf': '2'
}
rule feature_counts:
    input:
        expand(
            'results/{{dataset}}/star_align<star_method>/{sample_name}/Aligned.sortedByCoord.out.bam',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/star_align<star_method>/feature_counts.all.tsv'
    threads: 8
    params:
        samples = samples['sample_name'].tolist(),
        gtf = lambda wildcards: config[wildcards.dataset]['references']['gtf'],
        strandness = lambda wildcards: FEATURECOUNTS_STRANDNESS_MAP[config[wildcards.dataset]['lib']['strandness']]
    shell: 'featureCounts -p -s {params.strandness} --countReadPairs -B -T {threads} -a {params.gtf} -t exon -g gene_id -o {output} {input}'


rule deseq2_qc:
    input: rules.feature_counts.output
    output:
        cluster = 'results/{dataset}/deseq2/cluster<star_method>.pdf',
        pca = 'results/{dataset}/deseq2/pca<star_method>.pdf'
    params:
        samples = lambda wildcards: config[wildcards.dataset]['samples'],
        design = lambda wildcards: config[wildcards.dataset]['diff_qc']['design']
    script: '../scripts/deseq2_qc.R'


rule deseq2_diffexp:
    input:
        'results/{dataset}/star_align<star_method>/feature_counts.all.tsv'
    output:
        norm = 'results/{dataset}/deseq2/diffexp_norm<star_method>.{contrast}.{subset}.tsv',
        ma_plot = 'results/{dataset}/deseq2/diffexp_ma<star_method>.{contrast}.{subset}.pdf',
        volcano = 'results/{dataset}/deseq2/diffexp_volcano<star_method>.{contrast}.{subset}.pdf',
    params:
        samples = lambda wildcards: config[wildcards.dataset]['samples'],
        id_to_symbol = lambda wildcards: config[wildcards.dataset]['references']['id_to_symbol'],
        filter_by = lambda wildcards: config[wildcards.dataset]['diff_exp']['filter_by'][wildcards.subset],
        filter_value = '{subset}',
        contrast = lambda wildcards: config[wildcards.dataset]['diff_exp']['contrasts'][wildcards.contrast]
    script: '../scripts/deseq2_diffexp.R'
