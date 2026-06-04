#!/usr/bin/env python3

from snakemake.script import snakemake
import pandas as pd
import pyfaidx

GENOME_FASTA = pyfaidx.Fasta(snakemake.params.genome)

# U2: +GT/AG -CT/AC; U2 weak: +GC/AG -CT/GC; U12: +AT/AC -GT/AT;
# if strand info unkown, try to match one of the 6 motifs
# others are annotated as non-canonical
def check_splice_type(row):
    if row['strand'] == '+':
        if row['start_motif'] == 'GT' and row['end_motif'] == 'AG':
            return 'U2'
        elif row['start_motif'] == 'GC' and row['end_motif'] == 'AG':
            return 'U2w'
        elif row['start_motif'] == 'AT' and row['end_motif'] == 'AC':
            return 'U12'
        else:
            return 'NC'
    elif row['strand'] == '-':
        if row['start_motif'] == 'CT' and row['end_motif'] == 'AC':
            return 'U2'
        elif row['start_motif'] == 'CT' and row['end_motif'] == 'GC':
            return 'U2w'
        elif row['start_motif'] == 'GT' and row['end_motif'] == 'AT':
            return 'U12'
        else:
            return 'NC'
    else:
        motif = row['start_motif'] + '/' + row['end_motif']
        if motif in ['GT/AG', 'CT/AC']:
            return 'U2'
        elif motif in ['GC/AG', 'CT/GC']:
            return 'U2w'
        elif motif in ['AT/AC', 'GT/AT']:
            return 'U12'
        else:
            return 'NC'

# hisat2 splice site coordinate is 0-based and start at last exonic base;
# intronic region is [start+2, end]
# check intronic start and end bases.
splice_sites = pd.read_table(
    snakemake.input[0],
    header=None, names=['chr', 'start', 'end', 'strand'],
    dtype={'chr': str, 'start': int, 'end': int, 'strand': str}
)
# pyfaidx slicing uses 0-based coordinates, so intron splice motifs are at [start+1, start+3] and [end-2, end]
splice_sites['start_motif'] = splice_sites.apply(
    lambda x: GENOME_FASTA[x['chr']][x['start']+1:x['start']+3].seq.upper(), axis=1
)
splice_sites['end_motif'] = splice_sites.apply(
    lambda x: GENOME_FASTA[x['chr']][x['end']-2:x['end']].seq.upper(), axis=1
)

splice_sites['splice_type'] = splice_sites.apply(check_splice_type, axis=1)

if snakemake.params.coord_type == 'STAR':
    splice_sites['start'] += 2 # convert to STAR SJ.out coordinate

splice_sites.drop(columns=['start_motif', 'end_motif'], inplace=True)
splice_sites.to_csv(snakemake.output[0], sep='\t', index=False)
