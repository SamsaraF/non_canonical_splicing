import polars as pl
import numpy as np
import os
from pyfaidx import Fasta
from argparse import ArgumentParser
from typing import Literal


np.random.seed(int(os.getenv("NUMPY_RANDOM_SEED", 42)))
pl.set_random_seed(int(os.getenv("POLARS_RANDOM_SEED", 42)))


def extract_sequence(chr: str, strand: Literal['+','-'], exons: list[tuple[int, int]], genome_fasta: Fasta):
    """Return sequence of transcript. Exons are 1-based (from gtf)."""
    seq = ""
    for start, end in exons:
        start -= 1 # convert to 0-based for pyfaidx
        if strand == '+':
            seq += genome_fasta[chr][start:end].seq.upper()
        else:
            seq += genome_fasta[chr][start:end].reverse.complement.seq.upper()
    return seq


def format_fastq_read(read_id: str, seq: str):
    return f'@{read_id}\n{seq}\n+\n{"I" * len(seq)}'


def create_mock_fastq(sample_annot_path: str, fasta_path: str, usage_file_path: str, out_dir: str, read_length: int=100, error_rate: float=0.002):
    GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
    GENOME_FASTA = Fasta(fasta_path)
    sample_annot_lf = pl.scan_csv(
        sample_annot_path, separator='\t', comment_prefix='#', has_header=False,
        new_columns=GTF_COLUMNS
    )

    exons_lf = (
        sample_annot_lf
        .filter(pl.col('feature') == 'exon')
        .with_columns(
            transcript_id = pl.col('attr').str.extract(r'transcript_id "([^"]+)"'),
            coord = pl.concat_list('start', 'end')
        )
        .group_by('transcript_id', maintain_order=True)
        .agg(
            exons = pl.col('coord'),
            strand = pl.first('strand'),
            chr = pl.first('chr')
        )
        .with_columns(pl.col('exons').map_elements(lambda e: sorted(e, key=lambda x: x[0])))
        .select(
            ['transcript_id', 'chr', 'strand', 'exons'],
            sequence = pl.struct(
                'chr', 'strand', 'exons'
            ).map_elements(
                lambda row: extract_sequence(row['chr'], row['strand'], row['exons'], GENOME_FASTA),
                pl.String
            )
        )
        .with_columns(length = pl.col('sequence').str.len_chars())
    )
    transcript_exon_info_df = exons_lf.collect()

    transcript_usage_df = pl.read_csv(usage_file_path, separator='\t')
    sample_cols = [col for col in transcript_usage_df.columns if col.startswith('sample')]

    assert len(transcript_usage_df) == len(transcript_exon_info_df), "Usage file should have the same number of transcripts as annotation file"
    transcript_usage_df = (
        transcript_usage_df
        .join(transcript_exon_info_df, on='transcript_id', how='inner')
        # convert usage to read count
        .with_columns(
            **{sample: (pl.col(sample) * pl.col('length') / read_length).cast(pl.UInt32) for sample in sample_cols}
        )
    )

    for sample in sample_cols:
        out_fastq_path = f"{out_dir}/{sample}.fastq"
        with open(out_fastq_path, 'w') as sample_fastq:
            for row in transcript_usage_df.iter_rows(named=True):
                for read_id, start in enumerate(np.random.randint(0, row['length'] - read_length, row[sample])):
                    read_seq = row['sequence'][start: start + read_length]
                    # add random sequencing errors
                    error_mask = np.random.rand(len(read_seq)) < error_rate
                    random_seq_err = np.random.choice(list('ACGT'), size=len(read_seq))
                    read_seq_with_err = ''.join([random_seq_err[i] if error_mask[i] else read_seq[i] for i in range(len(read_seq))])

                    fastq_read = format_fastq_read(f"seed42:read_{read_id}", read_seq_with_err)
                    sample_fastq.write(fastq_read + '\n')


if __name__ == "__main__":
    parser = ArgumentParser(description="Create mock fastq files based on annotation and simulated transcript usage")
    parser.add_argument('-a', '--sample-annot', type=str, required=True, help='Path to sample annotation gtf file')
    parser.add_argument('-g', '--genome-path', type=str, required=True, help='Path to reference genome fasta file')
    parser.add_argument('-u', '--usage-file', type=str, required=True, help='Path to simulated transcript usage tsv file')
    parser.add_argument('-o', '--out-dir', type=str, required=True, help='Directory to output mock fastq files')
    parser.add_argument('--read_length', type=int, default=100, help='Read length for simulated fastq reads')
    parser.add_argument('--error_rate', type=float, default=0.002, help='Error rate for simulated sequencing errors')
    args = parser.parse_args()
    create_mock_fastq(args.sample_annot, args.genome_path, args.usage_file, args.out_dir, args.read_length, args.error_rate)
    