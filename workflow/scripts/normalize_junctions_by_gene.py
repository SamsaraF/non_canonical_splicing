#!/usr/bin/env python3
import pandas as pd
from snakemake.script import snakemake

junction_matrix = pd.read_table(
    snakemake.input['junction_matrix'], index_col=[0,1,2,3,4],
).astype(pd.Float64Dtype())

gene_matrix = pd.read_table(
    snakemake.input['gene_matrix'], index_col=0,
).astype(pd.Int32Dtype())

genes_bed = pd.read_table(
    snakemake.input['genes_bed'],
    index_col=3, usecols=[0,1,2], header=['chr', 'start', 'end']
)

def get_genes_by_junction(junc_chrod):
    chr, start, end, _ = junc_chrod
    genes_df = genes_bed[
        (genes_bed['chr'] == chr) &
        (genes_bed['start'] < start) &
        (genes_bed['end'] > end)
    ]
    gene_ids = genes_df.index.tolist()
    if len(gene_ids) == 0:
        gene_ids = None
    return gene_ids

junction_matrix['genes'] = junction_matrix.index.to_frame(index=False).apply(
    get_genes_by_junction, axis=1, result_type='expand'
)

junction_matrix.dropna(subset=['genes'], inplace=True)


def get_total_count_matrix(gene_ids):
    return gene_matrix.loc[gene_ids].sum()

total_count_matrix = junction_matrix['genes'].apply(get_total_count_matrix)

normalized_matrix = junction_matrix.divide(
    total_count_matrix.values, axis=0
).fillna(0)