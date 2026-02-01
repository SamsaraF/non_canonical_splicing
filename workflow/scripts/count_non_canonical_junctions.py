#!/usr/bin/env python3
import pandas as pd
from snakemake.script import snakemake


star_sjout_info_columns = {
    'chr': pd.StringDtype(),
    'start': pd.Int64Dtype(),
    'end': pd.Int64Dtype(),
    'strand': pd.Int64Dtype(), # 0 undefined, 1 +, 2 -
    'intron_motif': pd.Int64Dtype(), # 0 non-canonical, 1 GT/AG, 2 CT/AC, 3 GC/AG, 4 CT/GC, 5 AT/AC, 6 GT/AT
    'annotated': pd.Int64Dtype() # 0 no, 1 yes
}
star_sjout_data_columns = {
    'unique_reads': pd.Int64Dtype(),
    'multi_mapped_reads': pd.Int64Dtype(),
    'max_overhang': pd.Int64Dtype()
}
star_sjout_columns = star_sjout_info_columns | star_sjout_data_columns

sj_all = []
for f in snakemake.input:
    sample_sj = pd.read_table(
        f, index_col=[0,1,2,3,4,5],
        header=None, names=star_sjout_columns.keys(), dtype=star_sjout_columns
    )
    sample_non_canonical_count = sample_sj.xs(0, level='intron_motif')[['unique_reads']]
    sj_all.append(sample_non_canonical_count)

for sj_count, sample in zip(sj_all, snakemake.params.samples):
    sj_count.columns = [sample]

non_canonical_sj_matrix = pd.concat(sj_all, axis=1)
non_canonical_sj_matrix.fillna(0, inplace=True)

non_canonical_sj_matrix.to_csv(snakemake.output[0], sep='\t')
