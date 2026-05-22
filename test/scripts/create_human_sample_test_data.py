import os
import polars as pl


pl.set_random_seed(42)


def pick_from_human_genome(annot_path: str, gene_count=20, only_protein_coding=False):
    """Load gtf annotation and pick random genes to build a sample dataset."""
    annot_df = pl.scan_csv(
        annot_path, separator='\t', comment_prefix='#', has_header=False,
        new_columns=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
    )

    # count transcripts per gene, and filter genes with 3-20 transcripts
    transcript_query = (
        annot_df
        .select(['feature', 'attr'])
        .filter(pl.col('feature') == 'transcript')
        .filter(pl.col('attr').str.contains('gene_type "protein_coding";') if only_protein_coding else pl.lit(True))
        .with_columns(pl.col('attr').str.extract(r'gene_id "([^"]+)";').alias('gene_id'))
        .with_columns(
            pl.col('gene_id').unique(maintain_order=True),
            pl.col('gene_id').unique_counts().alias('transcript_count')
        )
        .filter(
            pl.col('transcript_count') >= 3,
            pl.col('transcript_count') <= 20
        )
    )
    sample_gene_pool = transcript_query.collect()

    # pick genes
    sample_gene_id = sample_gene_pool.sample(gene_count)['gene_id'].to_list()
    # get picked genes' annotation
    sample_query = (
        annot_df
        .filter(pl.col('attr').str.contains_any(sample_gene_id))
    )
    
    sample_df = sample_query.collect()
    
    # merge all genes into 1 chromosome
    sample_df

