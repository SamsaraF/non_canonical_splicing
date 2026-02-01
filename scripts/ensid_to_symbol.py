#!/usr/bin/env python3
# -*- coding=utf-8 -*-

import pandas as pd
import re


# bedops styled bed, 6 columns bed6 + 4 bedops extra columns
#genes_bed = 'references/mouse_ref/gencode.vM37.genes.bed'
#output_tsv = 'references/mouse_ref/ensid_to_symbol.tsv'
genes_bed = 'references/human_ref/gencode.v49.genes.bed'
output_tsv = 'references/human_ref/ensid_to_symbol.tsv'

genes_df = pd.read_table(
    genes_bed,
    header=None,
    names=['chr', 'start', 'end', 'id', 'score', 'strand', 'source', 'feature', 'frame', 'attr'],
)
genes_df.set_index('id', inplace=True, drop=False)

def extract_symbol(attr):
    return re.sub(r'.*gene_name "(\S*)";.*', r'\1', attr)

genes_df['symbol'] = genes_df['attr'].apply(extract_symbol)
symbol_df = genes_df[['id','symbol']].copy()

symbol_df.to_csv(output_tsv, sep='\t', index=False)
