import pandas as pd
import re


count_matrix = pd.read_csv(
    'results/hela_tm_tg/feature_counts_2pass_nc.all.tsv',
    sep='\t', comment='#', index_col=0, usecols=[0] + list(range(6, 18)),
)

count_matrix.rename(
    columns=lambda x: re.sub(
        r'results/hela_tm_tg/star_align_2pass_nc/(.*)/Aligned.sortedByCoord.out.bam$', 
        r'\1', x
    ),
    inplace=True,
)

count_matrix.to_csv('results/hela_tm_tg/feature_counts_2pass_nc.all.matrix.tsv', sep='\t')