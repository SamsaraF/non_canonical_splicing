#!/usr/bin/env python3

from sys import exit
from collections import defaultdict as dd
from argparse import ArgumentParser, FileType

# modified from hisat2_extract_splice_sites.py
# see: https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_extract_splice_sites.py


def extract_transcript_intronic_regions(gtf_file):
    genes = dd(list)
    trans = {}

    # Parse valid exon lines from the GTF file into a dict by transcript_id
    for line in gtf_file:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if '#' in line:
            line = line.split('#')[0].strip()

        try:
            chrom, source, feature, left, right, score, \
                strand, frame, values = line.split('\t')
        except ValueError:
            continue
        left, right = int(left), int(right)

        if feature != 'exon' or left >= right:
            continue

        values_dict = {}
        for attr in values.split(';'):
            if attr:
                attr, _, val = attr.strip().partition(' ')
                values_dict[attr] = val.strip('"')

        if 'gene_id' not in values_dict or \
                'transcript_id' not in values_dict:
            continue

        transcript_id = values_dict['transcript_id']
        if transcript_id not in trans:
            trans[transcript_id] = [chrom, strand, [[left, right]]]
            genes[values_dict['gene_id']].append(transcript_id)
        else:
            trans[transcript_id][2].append([left, right])

    # Sort exons and merge where separating introns are <=5 bps
    for transcript_id, [chrom, strand, exons] in trans.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[transcript_id] = [chrom, strand, tmp_exons, []]

        # calculate intronic regions
        for i in range(1, len(exons)):
            trans[transcript_id][3].append((exons[i-1][1] + 1, exons[i][0] - 1))
        
        # format intronic regions
        intron_regions = ';'.join(['{}-{}'.format(intron[0], intron[1]) for intron in trans[transcript_id][3]])
    
        print('{}\t{}\t{}\t{}'.format(transcript_id, chrom, strand, intron_regions))


if __name__ == '__main__':
    parser = ArgumentParser(
        description='Extract intronic regions for each transcript from a GTF file')
    parser.add_argument('gtf_file',
        nargs='?',
        type=FileType('r'),
        help='input GTF file (use "-" for stdin)')

    args = parser.parse_args()
    if not args.gtf_file:
        parser.print_help()
        exit(1)
    extract_transcript_intronic_regions(args.gtf_file)