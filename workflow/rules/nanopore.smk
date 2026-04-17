rule chopper:
    input: 'data/hela_tm_tg_ont/{sample_name}.fq.gz'
    output: 'results/hela_tm_tg_ont/fastq/{sample_name}.trim.fq.gz'
    threads: 4
    params:
        quality = 10,
        min_length = 300
    shell: 'gunzip -c {input} | chopper -q {params.quality} -l {params.min_length} -t {threads} | gzip > {output}'


rule quality_control:
    input: rules.chopper.output
    output: 'results/hela_tm_tg_ont/qc/{sample_name}/NanoPlot-report.html'
    threads: 4
    params:
        out_dir = subpath(output[0], parent=True)
    shell: 'NanoPlot -t {threads} --N50 --dpi 300 --fastq {input} --title {wildcards.sample_name} -o {params.out_dir}'


rule minimap2_index:
    params:
        genome = 'references/human_ref/GRCh38.p14.genome.fa',
        kmer = 14,
        window_size = 5
    output: 'references/human_ref/GRCh38.p14.genome.mmi'
    shell: 'minimap2 -k {params.kmer} -w {params.window_size} -d {output} {params.genome}'


rule minimap2:
    input: 
        fq = rules.chopper.output,
        idx = rules.minimap2_index.output
    output: temp('results/hela_tm_tg_ont/minimap2/{sample_name}.sam')
    threads: 4
    shell: 'minimap2 -a --splice -g2k -G200k -A1 -B2 -O2,24 -E1,0 -C0 -z200 -un --splice-flank=yes -t {threads} {input.idx} {input.fq} > {output}'


rule convert_to_bam:
    input: rules.minimap2.output
    output: temp('results/hela_tm_tg_ont/minimap2/{sample_name}.bam')
    shell: 'samtools view -bS {input} > {output}'


rule sort_bam:
    input: rules.convert_to_bam.output
    output: 'results/hela_tm_tg_ont/minimap2/{sample_name}.sorted.bam'
    threads: 4
    shell: 'samtools sort -@ {threads} -o {output} {input}'


rule index_minimap_bam:
    input: rules.sort_bam.output
    output: 'results/hela_tm_tg_ont/minimap2/{sample_name}.sorted.bam.bai'
    threads: 4
    shell: 'samtools index -@ {threads} {input}'


rule isoquant:
    input:
        bam = rules.sort_bam.output,
        idx = rules.index_minimap_bam.output
    output:
        model = 'results/hela_tm_tg_ont/isoquant/{sample_name}/{sample_name}.transcript_models.gtf',
        tpm = 'results/hela_tm_tg_ont/isoquant/{sample_name}/{sample_name}.discovered_transcript_tpm.tsv'
    params:
        genome = 'references/human_ref/GRCh38.p14.genome.fa',
        gtf = 'references/human_ref/gencode.v49.annotation.gtf',
        out_dir = 'results/hela_tm_tg_ont/isoquant/'
    shell: 'isoquant.py --reference {params.genome} --genedb {params.gtf} --complete_genedb --bam {input.bam} -d nanopore -o {params.out_dir} -p {wildcards.sample_name} --check_canonical --report_canonical all'


rule convert_novel_transcript_to_bed:
    input: rules.isoquant.output.model
    output: 'results/hela_tm_tg_ont/isoquant/{sample_name}/{sample_name}.transcript_novel.bed'
    shell: 'zsh workflow/scripts/isoquant_to_bed.sh -i {input} -o {output}'


rule assembly_ont_only:
    input:
        bam = rules.sort_bam.output,
        idx = rules.index_minimap_bam.output
    output: 'results/hela_tm_tg_ont/stringtie/{sample_name}.ONT_only.gtf'
    threads: 4
    params:
        gtf = 'references/human_ref/gencode.v49.annotation.gtf'
    shell: 'stringtie -L -p {threads} -G {params.gtf} -o {output} {input.bam}'


rule stringtie_ont_only_merge:
    input: expand('results/hela_tm_tg_ont/stringtie/{sample_name}.ONT_only.gtf', sample_name=samples['sample_name'].tolist())
    output: 'results/hela_tm_tg_ont/stringtie/merged.gtf'
    threads: 4
    params:
        gtf = 'references/human_ref/gencode.v49.annotation.gtf'
    shell: 'stringtie --merge -p {threads} -G {params.gtf} -o {output} {input}'


rule stringtie_ont_only_estimate:
    input:
        bam = rules.sort_bam.output,
        idx = rules.index_minimap_bam.output,
        merged_gtf = rules.stringtie_ont_only_merge.output
    output: 'results/hela_tm_tg_ont/stringtie/{sample_name}/estimate.gtf'
    threads: 4
    shell: 'stringtie -L -e -B -p {threads} -G {input.merged_gtf} -o {output} {input.bam}'
