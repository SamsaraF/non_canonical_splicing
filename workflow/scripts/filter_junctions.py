#!/usr/bin/env python3
import pandas as pd
from snakemake.script import snakemake


junction_matrix = pd.read_table(
    snakemake.input['junction_matrix'], index_col=[0,1,2,3,4],
).astype(pd.Float64Dtype())

if snakemake.params['filter_by'] == 'raw_count':
    junction_filtered = junction_matrix[
        junction_matrix
    ]