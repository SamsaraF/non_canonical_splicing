# -*- coding=utf-8 -*-

rule build_non_canonical_junction_matrix:
    input:
        expand(
            'results/{{dataset}}/star_align<star_method>/{sample_name}/SJ.out.tab',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/star_align<star_method>/junctions/non_canonical.count.tsv'
    params:
        samples = samples['sample_name'].tolist(),
        genes_bed = lambda wildcards: config[wildcards.dataset]['references']['genes_bed']
    script: '../scripts/count_non_canonical_junctions.py'


rule glob_total_mapped_counts:
    input:
        expand(
            'results/{{dataset}}/star_align<star_method>/{sample_name}/Log.final.out',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/star_align<star_method>/primary_mapped.stat.tsv'
    run:
        with open(output[0], 'w') as fout:
            fout.write('sample_name\ttotal_mapped_reads\n')
            for fpath in input:
                sample_name = fpath.split('/')[-2]
                with open(fpath) as fin:
                    for line in fin:
                        if 'Uniquely mapped reads number |' in line:
                            count = line.strip().split('\t')[1]
                            fout.write(f'{sample_name}\t{count}\n')
                            break

rule normalize_as_cpm:
    input:
        junction_matrix = rules.build_non_canonical_junction_matrix.output[0],
        mapped_counts = rules.glob_total_mapped_counts.output[0]
    output: 'results/{dataset}/star_align<star_method>/junctions/non_canonical.cpm.tsv'
    script: '../scripts/normalize_junctions_as_cpm.py'


rule junction_to_bed:
    input: rules.build_non_canonical_junction_matrix.output
    output: temp('results/{dataset}/star_align<star_method>/junctions/non_canonical.unsorted.bed')
    shell: "awk -v OFS='\\t' 'NR>1{{print $1, $2-1, $3}}' {input} > {output}"


rule sort_bed:
    input: rules.junction_to_bed.output
    output: 'results/{dataset}/star_align<star_method>/junctions/non_canonical.bed'
    shell: 'sort-bed {input} > {output}'


rule get_overlapped_genes:
    input: rules.sort_bed.output
    output: 'results/{dataset}/star_align<star_method>/junctions/non_canonical.overlapped_genes.tsv'
    params:
        genes_bed = 'references/human_ref/gencode.v49.genes.sorted.bed'
    shell: "bedmap --ec --delim '\\t' --echo --echo-map-id-uniq {input} {params.genes_bed} | awk -F '\\t' '$4' > {output}"


rule filter_junctions:
    input: rules.build_non_canonical_junction_matrix.output
    output: 'results/{dataset}/star_align<star_method>/junctions/non_canonical.count.filtered.tsv'
    params:
        filter_by = 'raw_count'
    script: '../scripts/filter_junctions.py'
