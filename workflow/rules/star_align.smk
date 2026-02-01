# -*- coding=utf-8 -*-

def get_read_command(wildcards, input):
    if input.r1.endswith('.gz'):
        return '--readFilesCommand zcat'
    else:
        return ''

rule star_align:
    input:
        unpack(get_trimmed_fastq)
    output:
        bam = 'results/{dataset}/<star_out>/{sample_name}/Aligned.out.bam',
        sort_bam = 'results/{dataset}/<star_out>/{sample_name}/Aligned.sortedByCoord.out.bam',
        gene_counts = 'results/{dataset}/<star_out>/{sample_name}/ReadsPerGene.out.tab',
        junctions = 'results/{dataset}/<star_out>/{sample_name}/SJ.out.tab',
        unmapped_r1 = 'results/{dataset}/<star_out>/{sample_name}/Unmapped.out.mate1',
        unmapped_r2 = 'results/{dataset}/<star_out>/{sample_name}/Unmapped.out.mate2'
    pathvars:
        star_out = 'star_align'
    params:
        max_intron = '1000000',
        ref = lambda wildcards: config[wildcards.dataset]['references']['star'],
        align_ends = 'Local',
        out_dir = subpath(output[0], parent=True),
        read_command = get_read_command
    threads: 16
    shell: 
        'STAR --runMode alignReads --outFilterType BySJout --outSAMattributes NH HI AS NM MD --runThreadN {threads} '
        '--outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --alignIntronMin 10 '
        '--outSJfilterOverhangMin 2 2 2 2 --outSJfilterDistToOtherSJmin 0 0 5 10 --outSJfilterCountTotalMin 1 1 1 1 '
        '--alignIntronMax {params.max_intron} --alignMatesGapMax {params.max_intron} ' # mannual set max intron length
        '--alignEndsType {params.align_ends} '
        #'--outSAMattrIHstart 0 --outSAMstrandField intronMotif ' # cufflinks compatibility
        '--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 '
        ' --outSAMtype BAM Unsorted SortedByCoordinate '
        '--scoreGapNoncan -4 '
        '--quantMode TranscriptomeSAM GeneCounts --outReadsUnmapped Fastx --genomeDir {params.ref} '
        '--readFilesIn {input.r1} {input.r2} {params.read_command} --outFileNamePrefix {params.out_dir}/'


use rule star_align as star_align_soft_clip with:
    input:
        r1 = 'results/{dataset}/star_align/{sample_name}/Unmapped.out.mate1',
        r2 = 'results/{dataset}/star_align/{sample_name}/Unmapped.out.mate2'
    pathvars:
        star_out = 'star_align_soft_clip'
