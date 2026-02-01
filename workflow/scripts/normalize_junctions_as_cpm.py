#!/usr/bin/env python3
import pandas as pd
from snakemake.script import snakemake

junction_matrix = pd.read_table(
    snakemake.input['junction_matrix'], index_col=[0,1,2,3,4],
).astype(pd.Float64Dtype())

total_primary_counts = pd.read_table(
    snakemake.input['mapped_counts'], index_col=0,
)
norm_factors = total_primary_counts['total_mapped_reads'].to_dict()

cpm_matrix = (junction_matrix * 1e6).div(
    norm_factors, axis=1
)

cpm_matrix.to_csv(snakemake.output[0], sep='\t')
