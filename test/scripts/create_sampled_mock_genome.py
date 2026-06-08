import polars as pl
from typing import Literal
from pyfaidx import Fasta
from argparse import ArgumentParser
import os


pl.set_random_seed(int(os.getenv("POLARS_RANDOM_SEED", 42)))


def calc_gap_length(gene_df: pl.DataFrame, current_pos: int, chr_size: int, compared: Literal['next', 'prev']='next', max_length=100000) -> int:
    """Tool function, limit max gap length to 100kb."""
    if compared not in ['next', 'prev']:
        raise ValueError("compared must be 'next' or 'prev'")
    if gene_df.is_empty():
        print(f"Warning: no gene found.")
        if compared == 'next':
            return min(chr_size - current_pos, max_length)
        else:
            return min(current_pos - 1, max_length)
    elif gene_df.shape[0] == 1:
        if compared == 'next':
            return min(gene_df[0, 'start'] - current_pos - 1, max_length)
        else:
            return min(current_pos - gene_df[0, 'end'] - 1, max_length)
    else:
        raise ValueError("gene_df should contain one gene")
    

def write_fasta(out_file_path: str, sequences: dict[str, str], width: int=80):
    """Write sequences to FASTA file with specified width."""
    with open(out_file_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>{seq_id}\n')
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + '\n')


def filter_attributes(attr: str) -> str:
    """Filter GTF attribute string to keep only specified keys."""
    KEPT_KEYS = ['gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_type', 'level']
    attr_dict = {}
    for item in attr.split(';'):
        item = item.strip()
        if not item:
            continue
        key, value = item.split(' ', 1)
        if key in KEPT_KEYS:
            attr_dict[key] = value
    return '; '.join(f'{key} {value}' for key, value in attr_dict.items()) + ';'


def create_sampled_mock_genome(
    annot_path: str, genome_fasta_path: str, out_path,
    gene_count=20, only_protein_coding=False,
    transcript_count_min=4, transcript_count_max=20, gene_length_min=4000, gene_length_max=40000,
    max_intergenic_length=100000
):
    """
    Load gtf annotation and pick random genes to build a sample dataset.
    All picked genes are merged to 1 mock chromosome, connected by itself tailing intergenic region.
    Output a gtf file and a fasta file for the reconstructed chromosome.
    """

    GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
    annot_lf = pl.scan_csv(
        annot_path, separator='\t', comment_prefix='#', has_header=False,
        new_columns=GTF_COLUMNS
    )

    # count transcripts per gene, and filter genes with too many or too few transcripts
    transcript_lf = (
        annot_lf
        .select(['feature', 'attr'])
        .filter(pl.col('feature') == 'transcript')
        .filter(pl.col('attr').str.contains('gene_type "protein_coding";') if only_protein_coding else pl.lit(True))
        .with_columns(gene_id = pl.col('attr').str.extract(r'gene_id "([^"]+)";'))
        .select(
            pl.col('gene_id').unique(maintain_order=True),
            pl.col('gene_id').unique_counts().alias('transcript_count')
        )
        .filter(pl.col('transcript_count').is_between(transcript_count_min, transcript_count_max))
    )
    filtered_gene_pool = transcript_lf.collect()
    filtered_gene_ids = filtered_gene_pool['gene_id'].to_list()

    # calc gene length, and filter genes with too long or too short length
    gene_lf = (
        annot_lf
        .filter(pl.col('feature') == 'gene')
        .with_columns(gene_id = pl.col('attr').str.extract(r'gene_id "([^"]+)";'))
        .filter(pl.col('gene_id').is_in(filtered_gene_ids))
        .with_columns(gene_length = pl.col('end') - pl.col('start') + 1)
        .filter(pl.col('gene_length').is_between(gene_length_min, gene_length_max))
    )

    # pick random genes
    sample_gene_ids = gene_lf.collect().sample(gene_count)['gene_id'].to_list()
    sample_lf = (
        annot_lf
        .with_columns(gene_id = pl.col('attr').str.extract(r'gene_id "([^"]+)";'))
        .filter(pl.col('gene_id').is_in(sample_gene_ids))
        .select(GTF_COLUMNS + ['gene_id'])
    )

    sample_annot_df = sample_lf.collect()
    sample_gene_df = (
        sample_annot_df
        .filter(pl.col('feature') == 'gene')
        .with_columns(gene_length = pl.col('end') - pl.col('start') + 1)
    )

    sample_chrs = sample_gene_df['chr'].unique().to_list()
    all_needed_genes_lf = (
        annot_lf
        .filter(pl.col('feature') == 'gene', pl.col('chr').is_in(sample_chrs))
        .select(pl.all().sort_by('start').over('chr'))
    )
    all_needed_genes_df = all_needed_genes_lf.collect()

    total_length = 0
    coord_shift_records = {} # genome coordinate start (bp, 1-based) for each gene after shifting
    gene_id_mapping = {} # mapping from original gene id to new gene id (gene_1, gene_2, ...)
    genome_fasta = Fasta(genome_fasta_path) # slicing use 0-based
    mock_genome_sequence = ''

    for idx, row in enumerate(sample_gene_df.iter_rows(named=True)):
        chr_size = len(genome_fasta[row['chr']])
        gene_id_mapping[row['gene_id']] = idx + 1
        if idx == 0:
            prev_gene_df = (
                all_needed_genes_df
                .filter(pl.col('chr') == row['chr'], pl.col('end') < row['start'])
                .tail(1)
            )
            prev_gap_length = calc_gap_length(prev_gene_df, row['start'], chr_size, 'prev', max_intergenic_length)
            mock_genome_sequence += genome_fasta[row['chr']][(row['start'] - prev_gap_length - 1):(row['start'] - 1)].seq.upper()
            total_length += prev_gap_length
        
        coord_shift_records[row['gene_id']] = row['start'] - total_length - 1
        # get next gene coord to calculate intergenic region length
        next_gene_df = (
            all_needed_genes_df
            .filter(pl.col('chr') == row['chr'], pl.col('start') > row['end'])
            .head(1)
        )
        next_gap_length = calc_gap_length(next_gene_df, row['end'], chr_size, 'next', max_intergenic_length)
        mock_genome_sequence += genome_fasta[row['chr']][(row['start'] - 1):(row['end'] + next_gap_length)].seq.upper()
        total_length += (row['gene_length'] + next_gap_length)

    mock_annot_df = (
        sample_annot_df
        .with_columns(
            chr = pl.lit('chr1'),
            source = pl.lit('TEST'),
            gene_index = pl.col('gene_id').replace_strict(gene_id_mapping),
            shift_value = pl.col('gene_id').replace_strict(coord_shift_records)
        )
        .with_columns(
            start = pl.col('start') - pl.col('shift_value'),
            end = pl.col('end') - pl.col('shift_value'),
            # execute gene_id replacement
            attr = pl.col('attr').str.replace(r'gene_id "([^"]+)";', pl.format('gene_id "gene_{}";', pl.col('gene_index')))
        )
    )

    # extract transcript id and rename it to trans_<gene_id>.1, trans_<gene_id>.2, ...
    transcript_id_df = (
        mock_annot_df
        .filter(pl.col('feature') == 'transcript')
        .with_columns(transcript_id = pl.col('attr').str.extract(r'transcript_id "([^"]+)";'))
        .select(pl.col('transcript_id'), pl.col('gene_index'))
        .unique(maintain_order=True)
        .with_columns(transcript_index = pl.col('transcript_id').cum_count().over('gene_index'))
        .with_columns(new_transcript_id = pl.format('transcript_{}.{}', pl.col('gene_index'), pl.col('transcript_index')))
        .select('transcript_id', 'new_transcript_id')
    )
    transcript_id_mapping = dict(transcript_id_df.iter_rows())

    mock_annot_df = (
        mock_annot_df
        # execute replacement
        .with_columns(
            pl.col('attr').str.replace_many(list(transcript_id_mapping.keys()), list(transcript_id_mapping.values())),
        )
        # keep only gtf columns for output
        .select(GTF_COLUMNS)
        # remove optional attribute fields
        .with_columns(pl.col('attr').map_elements(filter_attributes))
    )

    mock_annot_df.write_csv(f'{out_path}/annotation.gtf', include_header=False, separator='\t', quote_style='never')
    write_fasta(f'{out_path}/genome.fa', {'chr1': mock_genome_sequence})


if __name__ == '__main__':
    parser = ArgumentParser(description='Create mock genome for testing by random sampling genes from provided genome.')
    parser.add_argument('-a', '--annot', type=str, required=True, help='Path to the GENCODE style annotation GTF file.')
    parser.add_argument('-g', '--genome', type=str, required=True, help='Path to the genome FASTA file.')
    parser.add_argument('-o', '--out-dir', type=str, required=True, help='Output directory for the generated GTF and FASTA files.')
    parser.add_argument('--gene-count', type=int, default=20, help='Number of genes to sample for the mock genome.')
    parser.add_argument('--only-protein-coding', action='store_true', help='Whether to only include protein-coding genes in the sampling.')

    args = parser.parse_args()

    if os.path.exists(f'{args.out_dir}/annotation.gtf') or os.path.exists(f'{args.out_dir}/genome.fa'):
        option = input(f"Output files already exist in {args.out_dir}. Do you want to overwrite them? (y/n): ")
        if option.lower() != 'y':
            print("Aborting.")
            exit(0)
    
    create_sampled_mock_genome(
        annot_path=args.annot,
        genome_fasta_path=args.genome,
        out_path=args.out_dir,
        gene_count=args.gene_count,
        only_protein_coding=args.only_protein_coding,
        transcript_count_min=4,
        transcript_count_max=20,
        gene_length_min=4000,
        gene_length_max=40000,
        max_intergenic_length=100000
    )
