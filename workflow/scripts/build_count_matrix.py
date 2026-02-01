#!/usr/bin/env python3
import pandas as pd
from snakemake.script import snakemake


def get_star_readscount_column(strandness):
    if strandness == 'un':
        return 1
    elif strandness == 'rf':
        return 2
    elif strandness == 'fr':
        return 3
    else:
        raise ValueError(f"Unknown strandness: {strandness}. Supported values are 'un', 'rf', 'fr'.")
    

strandness = snakemake.params.get('strandness', 'un')
counts = [
    pd.read_table(f, index_col=0, usecols=[0, get_star_readscount_column(strandness)], header=None, skiprows=4)
    for f in snakemake.input
]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

count_matrix = pd.concat(counts, axis=1)
count_matrix.index.name = 'gene_id'

count_matrix.to_csv(snakemake.output[0], sep='\t')
