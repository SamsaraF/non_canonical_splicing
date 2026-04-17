# -*- coding=utf-8 -*-

configfile: 'config/config.yaml'

import pandas as pd

# Load sample information
samples = pd.read_csv(
    config[DATASET]['samples'], sep='\t', dtype={'sample_name': str}
).set_index('sample_name', drop=False).sort_index()
samples['replicate'] = samples.groupby(['condition_1', 'condition_2'])['sample_name'].cumcount() + 1
samples_reverse = samples.set_index(['condition_1', 'condition_2', 'replicate'], drop=False).sort_index()


def get_raw_fastq(wildcards):
    if pd.isna(samples.loc[wildcards.sample_name, 'fq1']):
        return {
            'r1': f'results/{wildcards.dataset}/fastq/{wildcards.sample_name}_1.fastq.gz',
            'r2': f'results/{wildcards.dataset}/fastq/{wildcards.sample_name}_2.fastq.gz'
        }
    else:
        fq1_path = samples.loc[wildcards.sample_name, 'fq1']
        fq2_path = samples.loc[wildcards.sample_name, 'fq2']
        return {
            'r1': f'data/{wildcards.dataset}/{fq1_path}',
            'r2': f'data/{wildcards.dataset}/{fq2_path}'
        }


def get_trimmed_fastq(wildcards):
    if config[wildcards.dataset]['lib']['adapter'] != 'none':
        return {
            'r1': f'results/{wildcards.dataset}/cleaned_fastq/{wildcards.sample_name}_val_1.fq.gz',
            'r2': f'results/{wildcards.dataset}/cleaned_fastq/{wildcards.sample_name}_val_2.fq.gz'
        }
    else:
        return get_raw_fastq(wildcards)


def get_sorted_bam(wildcards):
    query = (wildcards.condition_1, wildcards.condition_2, int(wildcards.replicate))
    sample_name = samples_reverse.loc[query, 'sample_name']
    return f'results/{wildcards.dataset}/star_align<star_method>/{sample_name}/Aligned.sortedByCoord.out.bam'
