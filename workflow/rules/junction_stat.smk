# -*- coding=utf-8 -*-

rule build_non_canonical_junction_matrix:
    input:
        expand(
            'results/{{dataset}}/star_align/{sample_name}/SJ.out.tab',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/junctions/non_canonical_junctions.all.tsv'
    params:
        samples = samples['sample_name'].tolist()
    script: '../scripts/count_non_canonical_junctions.py'


rule glob_total_mapped_counts:
    input:
        expand(
            'results/{{dataset}}/star_align/{sample_name}/Log.final.out',
            sample_name=samples['sample_name'].tolist()
        )
    output: 'results/{dataset}/star_align/primary_mapped.stat.tsv'
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
        junction_matrix = 'results/{dataset}/junctions/non_canonical_junctions.all.tsv',
        mapped_counts = 'results/{dataset}/star_align/primary_mapped.stat.tsv'
    output: 'results/{dataset}/junctions/non_canonical_junctions.all.cpms.tsv'
    script: '../scripts/normalize_junctions_as_cpm.py'