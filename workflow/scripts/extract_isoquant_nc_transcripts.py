import pandas as pd
import re
import matplotlib.pyplot as plt


def is_non_canonical(row: pd.Series) -> bool:
    if re.search('Canonical "False";', row['attribute']):
        return True
    return False


def load_non_canonical_transcripts(gtf_path: str) -> pd.DataFrame:
    """
    Load non-canonical transcripts from a GTF file into a pandas DataFrame.
    kept columns: ID(col 1), source(col 2), feature(col 3), start(col 4), end(col 5), strand(col 7), attribute(col 9)

    Args:
        gtf_path (str): The path to the GTF file.

    Returns:
        pd.DataFrame: The loaded DataFrame.
    """
    gtf = pd.read_csv(
        gtf_path, sep='\t', comment='#', header=None, usecols=[0, 1, 2, 3, 4, 6, 8],
        names=['chromosome', 'source', 'feature', 'start', 'end', 'strand', 'attribute'],
    )
    transcripts = gtf[gtf['feature'] == 'transcript'].copy()
    nc_transcripts = transcripts[transcripts.apply(is_non_canonical, axis=1)].copy()
    nc_transcripts.reset_index(inplace=True)
    return nc_transcripts


def load_novel_transcripts_bed(bed_path: str) -> pd.DataFrame:
    bed = pd.read_csv(
        bed_path, sep='\t', header=None, usecols=[0, 1, 2, 3, 5, 9, 10, 11],
        names=['chromosome', 'start', 'end', 'name', 'strand', 'block_count', 'block_sizes', 'block_starts'],
    )
    return bed


DATA_DIR = 'results/hela_tm_tg_ont/isoquant/'
KO_SAMPLE = 'KO_Tg'
WT_SAMPLE = 'WT_Tg'

ko_novel = load_novel_transcripts_bed(f'{DATA_DIR}/{KO_SAMPLE}/{KO_SAMPLE}.transcript_novel.bed')
wt_novel = load_novel_transcripts_bed(f'{DATA_DIR}/{WT_SAMPLE}/{WT_SAMPLE}.transcript_novel.bed')
all_novel = pd.merge(
    ko_novel, wt_novel,
    on=['chromosome', 'start', 'end', 'strand', 'block_count', 'block_sizes', 'block_starts'],
    how='outer', suffixes=('_ko', '_wt'), indicator='source_info'
)
all_novel['name'] = all_novel.apply(lambda row: f"transcript.{row.index}", axis=1)

wt_name_map = {}
ko_name_map = {}
for idx, row in all_novel.iterrows():
    if row['source_info'] == 'left_only':
        wt_name_map[row['name_wt']] = row['name']
    elif row['source_info'] == 'right_only':
        ko_name_map[row['name_ko']] = row['name']
    else:
        wt_name_map[row['name_wt']] = row['name']
        ko_name_map[row['name_ko']] = row['name']


ko_tpm = pd.read_csv(f'{DATA_DIR}/{KO_SAMPLE}/{KO_SAMPLE}.discovered_transcript_tpm.tsv', sep='\t')
wt_tpm = pd.read_csv(f'{DATA_DIR}/{WT_SAMPLE}/{WT_SAMPLE}.discovered_transcript_tpm.tsv', sep='\t')

ko_tpm['new_name'] = ko_tpm['#feature_id'].map(ko_name_map)
ko_tpm['new_name'].fillna(ko_tpm['#feature_id'], inplace=True)
wt_tpm['new_name'] = wt_tpm['#feature_id'].map(wt_name_map)
wt_tpm['new_name'].fillna(wt_tpm['#feature_id'], inplace=True)

# build plot data
plot_data = pd.merge(
    ko_tpm[['new_name', 'TPM']], wt_tpm[['new_name', 'TPM']],
    on='new_name', how='outer', suffixes=('_ko', '_wt')
)

