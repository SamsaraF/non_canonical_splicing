# -*- coding=utf-8 -*-


rule stringtie_mix_assembly:
    input:
        nano_bam = 'results/hela_tm_tg_ont/minimap2/{sample_name}.sorted.bam',
        rna_bam = 'results/hela_tm_tg/star_align_2pass_nc/{sample_name}_1/Aligned.sortedByCoord.out.bam'
    output: 'results/hela_tm_tg_mix/stringtie/{sample_name}/1.assemble.gtf'
    params:
        gff = 'references/human_ref/gencode.v49.annotation.gff3',
        use_viral = ''
    threads: 4
    shell: 'stringtie {params.use_viral} --mix -p {threads} -G {params.gff} -o {output} {input.rna_bam} {input.nano_bam}'


use rule stringtie_mix_assembly as stringtie_mix_viral_assembly with:
    output: 'results/hela_tm_tg_mix/stringtie_viral/{sample_name}/1.assemble.gtf'
    params:
        use_viral = '--viral'


rule stringtie_mix_viral_merge:
    input: expand('results/hela_tm_tg_mix/stringtie_viral/{sample_name}/1.assemble.gtf', sample_name=samples['sample_name'].tolist())
    output: 'results/hela_tm_tg_mix/stringtie_viral/merged.gtf'
    params:
        gff = 'references/human_ref/gencode.v49.annotation.gff3'
    threads: 4
    shell: 'stringtie --merge -p {threads} -G {params.gff} -o {output} {input}'


rule extarct_mix_viral_splice_sites:
    input: rules.stringtie_mix_viral_merge.output
    output: 'results/hela_tm_tg_mix/stringtie_viral/merged.splice_sites.tsv' # 0-based, start at last exonic base
    shell: 'hisat2_extract_splice_sites.py {input} > {output}'


rule append_mix_viral_ss_type:
    input: rules.extarct_mix_viral_splice_sites.output
    output: 'results/hela_tm_tg_mix/stringtie_viral/merged.splice_sites.with_types.tsv'
    params:
        genome = 'references/human_ref/GRCh38.p14.genome.fa',
        coord_type ='STAR' # to STAR junction coord, 1-based, start at first intronic base
    script: '../scripts/check_splice_sites_type.py'

