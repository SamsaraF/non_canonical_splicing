from argparse import ArgumentParser
import polars as pl
import numpy as np


pl.set_random_seed(42)
np.random.seed(42)


SPLICE_LENGTH_DISTRIBUTION = [16, 26, 32, 64, 128]


def find_end(exons: list[dict[str, int]], start: int, length: int) -> int:
    need = length
    cur = start
    for exon in exons:
        if exon['end'] < cur:
            continue
        if cur < exon['start']:
            cur = exon['start']
        available = exon['end'] - cur + 1
        if available >= need:
            return int(cur + need - 1)
        need -= available
        cur = exon['end'] + 1
    raise ValueError("Not enough length in intervals")


def remove_splice_site(exons: list[dict[str, int]], splice_start: int, splice_end: int) -> list[dict[str, int]]:
    result = []
    for exon in exons:
        if exon['end'] < splice_start:
            result.append(exon)
        elif exon['start'] > splice_end:
            result.append(exon)
        else:
            if exon['start'] <= splice_start - 1:
                result.append({'start': exon['start'], 'end': splice_start - 1, 'length': splice_start - exon['start']})
            if exon['end'] >= splice_end + 1:
                result.append({'start': splice_end + 1, 'end': exon['end'], 'length': exon['end'] - splice_end})
    return result


def generate_random_splice(exons: list[dict[str, int]], length_distr: list[int]=SPLICE_LENGTH_DISTRIBUTION) -> tuple[int, int]:
    """
    Generate random splice start and end by input exon regions.
    - exons: list of exon records sorted by start position, like [{'start': 100, 'end': 200, 'length': 101}, ...]
    """
    splice_length: int = np.random.choice(length_distr)
    total_length = sum(exon['length'] for exon in exons)
    # randomly pick splice start pos
    relative_start = np.random.randint(1, total_length - splice_length + 2) # add 1 to upperbound to make inclusive
    
    splice_start = find_end(exons, exons[0]['start'], relative_start)
    splice_end = find_end(exons, splice_start, splice_length)

    return (splice_start, splice_end)


def add_splicing_site(annot_path: str, out_path: str, event_rate: float=0.5):
    """
    Pick random gene from annotation, and add random deletions to 1 transcript to mock novel splicing sites.
    """

    GTF_COLUMNS = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr']
    annot_lf = pl.scan_csv(
        annot_path, separator='\t', comment_prefix='#', has_header=False,
        new_columns=GTF_COLUMNS
    )
    sample_genes_lf = (
        annot_lf
        .filter(pl.col('feature') == 'gene')
        .select(gene_id = pl.col('attr').str.extract(r'gene_id "([^"]+)"'))
        .select('gene_id')
    )
    all_genes = sample_genes_lf.collect()['gene_id'].to_list()
    sample_genes = np.random.choice(all_genes, size=int(len(all_genes) * event_rate), replace=False)

    novel_annotations = []
    for gene_id in all_genes:
        gene_elements_lf = annot_lf.filter(pl.col('attr').str.contains(f'gene_id "{gene_id}"'))
        old_elements_df = gene_elements_lf.collect()

        if gene_id not in sample_genes:
            novel_annotations.append(old_elements_df)
        else:
            # pick 1 transcript randomly and extract exon region
            all_transcript_ids_lf = (
                gene_elements_lf
                .filter(pl.col('feature') == 'transcript')
                .select(transcript_id = pl.col('attr').str.extract(r'transcript_id "([^"]+)"'))
            )
            random_transcript = all_transcript_ids_lf.collect().sample(n=1)
            picked_transcript_id = random_transcript['transcript_id'][0]
            picked_transcript_elements_lf = gene_elements_lf.filter(pl.col('attr').str.contains(f'transcript_id "{picked_transcript_id}"'))

            exons_lf = (
                picked_transcript_elements_lf
                .filter(pl.col('feature') == 'exon')
                .with_columns(length = pl.col('end') - pl.col('start') + 1)
                .select('start', 'end', 'length')
                .sort('start')
            )
            exon_regions = exons_lf.collect().to_dicts()
            splice_start, splice_end = generate_random_splice(exon_regions)

            new_exon_regions = remove_splice_site(exon_regions, splice_start, splice_end)

            transcript_record = (
                picked_transcript_elements_lf
                .filter(pl.col('feature') == 'transcript')
                .with_columns(attr = pl.col('attr').str.replace(f'transcript_id "{picked_transcript_id}"', f'transcript_id "{picked_transcript_id}.N"'))
            ).collect()
            novel_transcript_df = transcript_record.clone()
            assert len(transcript_record) == 1
            for region in new_exon_regions:
                region_df = (
                    transcript_record
                    .with_columns(
                        feature = pl.lit('exon'),
                        start = pl.lit(region['start']).cast(pl.Int64),
                        end = pl.lit(region['end']).cast(pl.Int64)
                    )
                )
                novel_transcript_df = novel_transcript_df.vstack(region_df)

            # # re-calculate regions
            # # handle splice site cross multiple exons first
            # novel_transcript_lf = (
            #     picked_transcript_elements_lf
            #     # remove elements fully covered by splice site
            #     .filter(~((splice_start < pl.col('start')) & (splice_end > pl.col('end'))))
            #     # all possible splice are not 3*n, there must be frameshift
            #     # remove all CDS/start_codon/stop_codon features to avoid complexity of re-assigning them
            #     .filter(~pl.col('feature').is_in(['CDS', 'start_codon', 'stop_codon']))
            #     # change start/end of elements partially covered by splice site
            #     .with_columns(
            #         start = pl.when((pl.col('start') < splice_end) & (splice_end < pl.col('end')) & (splice_start < pl.col('start')))
            #         .then(pl.lit(splice_end + 1))
            #         .otherwise(pl.col('start')),
            #         end = pl.when((pl.col('start') < splice_start) & (splice_start < pl.col('end')) & (splice_end > pl.col('end')))
            #         .then(pl.lit(splice_start - 1))
            #         .otherwise(pl.col('end')),
            #         # rename transcript id
            #         attr = pl.col('attr').str.replace(f'transcript_id "{picked_transcript_id}"', f'transcript_id "{picked_transcript_id}.N"')
            #     )
            # )
            # # consider splice site only on 1 exon
            # fully_covered_exon_lf = (
            #     novel_transcript_lf
            #     .filter((pl.col('start') < splice_start) & (splice_end < pl.col('end')))
            # )

            novel_annotations.append(pl.concat([old_elements_df, novel_transcript_df]))

    output_df: pl.DataFrame = pl.concat(novel_annotations)
    output_df.write_csv(out_path, include_header=False, separator='\t', quote_style='never')


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('-a', '--annot-path', required=True, help='Path to input annotation gtf file')
    parser.add_argument('-o', '--out-path', required=True, help='Path to output gtf file with novel splicing sites')
    parser.add_argument('-e', '--event-rate', type=float, default=0.5, help='Fraction of genes to modify')
    args = parser.parse_args()

    add_splicing_site(args.annot_path, args.out_path, args.event_rate)
