import polars as pl
import numpy as np
from argparse import ArgumentParser


np.random.seed(42)
pl.set_random_seed(42)


def simulate_transcript_abundances(
    splice_annot_path: str, out_path: str, replicate_count: int = 5
):
    GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
    annot_lf = pl.scan_csv(
        splice_annot_path, separator='\t', comment_prefix='#', has_header=False,
        new_columns=GTF_COLUMNS
    )
    transcript_count_lf = (
        annot_lf
        .filter(pl.col('feature') == 'transcript')
        .select(
            gene_id=pl.col('attr').str.extract(r'gene_id "([^"]+)"'),
            transcript_id=pl.col('attr').str.extract(r'transcript_id "([^"]+)"')
        )
        .unique(subset='transcript_id', maintain_order=True)
        .group_by('gene_id', maintain_order=True)
        .agg(
            transcript_count=pl.len(),
            transcript_ids=pl.col('transcript_id')
        )
    )

    transcript_count_df = transcript_count_lf.collect()
    abundance_df = transcript_count_df.with_columns(
        gene_abundance = pl.lit((np.random.normal(100, 20, size=transcript_count_df.height))).cast(pl.UInt32)
    ) # simulate gene total abundances as base value

    sample_labels = [f'sample.control_{i}' for i in range(1, replicate_count + 1)] + [f'sample.experiment_{i}' for i in range(1, replicate_count + 1)]
    control_libsize = np.random.normal(10, 1.5, replicate_count)
    experiment_libsize = np.random.normal(10, 1.5, replicate_count)

    out_rows_list = []
    for row in abundance_df.iter_rows(named=True):
        transcript_ids = row['transcript_ids']
        # half for experiment, half for control. assume transcripts are nearly equally expressed
        control_weights = [1 if idx.endswith('.N') else 20 for idx in transcript_ids ]
        experiment_weights = [20] * len(transcript_ids) # assume novel splicing over-expressed in experiment

        transcript_usages = np.random.dirichlet(control_weights + experiment_weights) * row['gene_abundance'] * 2

        control_matrix = np.outer(transcript_usages[:len(transcript_ids)], control_libsize)
        experiment_matrix = np.outer(transcript_usages[len(transcript_ids):], experiment_libsize)

        all_matrix = np.hstack([control_matrix, experiment_matrix])
        # add random noise
        all_matrix *= np.random.normal(1, 0.1, all_matrix.shape)

        for i, transcript_id in enumerate(transcript_ids):
            out_rows_list.append({
                'gene_id': row['gene_id'],
                'transcript_id': transcript_id,
                **{sample_labels[j]: int(all_matrix[i, j]) for j in range(all_matrix.shape[1])}
            })

    out_df = pl.DataFrame(out_rows_list, schema={
        'gene_id': pl.Utf8,
        'transcript_id': pl.Utf8,
        **{label: pl.UInt32 for label in sample_labels}
    })
    out_df.write_csv(out_path, separator='\t')


if __name__ == '__main__':
    parser = ArgumentParser(description="Simulate transcript abundances for control and experiment groups based on a splicing annotation.")
    parser.add_argument('-a', '--splice_annot_path', type=str, required=True, help='Path to the splicing annotation GTF file.')
    parser.add_argument('-o', '--out_path', type=str, required=True, help='Path to output the simulated transcript abundances TSV file.')
    parser.add_argument('-n', '--replicate_count', type=int, default=5, help='Number of replicates for each group (control and experiment).')
    
    args = parser.parse_args()
    simulate_transcript_abundances(args.splice_annot_path, args.out_path, args.replicate_count)
