#!/usr/bin/env python3
# -*- coding=utf-8 -*-


from argparse import ArgumentParser
from collections import defaultdict
import pysam


def polish_event():
    pass


DEL_THRESHOLD = 10
OVERHANG_THRESHOLD = 1


def extract_splice_events(bam_file: pysam.AlignmentFile, del_threshold=DEL_THRESHOLD, overhang_threshold=OVERHANG_THRESHOLD):
    splicing_event_count = defaultdict(int)

    for read in bam_file.fetch():
        if read.is_unmapped:
            continue
        if read.cigarstring and 'N' in read.cigarstring:
            chrom = bam_file.get_reference_name(read.reference_id)
            start = read.reference_start
            for cigar in read.cigartuples:
                if cigar[0] in [0, 2, 3, 7, 8]: # M, D, N, =, X
                    # detect splice or del > DEL_THRESHOLD
                    if cigar[0] == 3 or (cigar[0] == 2 and cigar[1] >= del_threshold):
                        splicing_event_count[(chrom, start, start + cigar[1])] += 1
                    start += cigar[1]

    return splicing_event_count


if __name__ == '__main__':
    parser = ArgumentParser(description='Detect splice events from BAM file')
    parser.add_argument('-b', '--bam', required=True, help='Path to the BAM file')
    parser.add_argument('-o', '--output', required=True, help='Path to save the output splice events')
    args = parser.parse_args()

    
    with open(args.output, 'w') as out_f:
        bamfile = pysam.AlignmentFile(args.bam, "rb")
        splicing_event_count = extract_splice_events(bamfile)
        bamfile.close()

        splicing_count_sorted = sorted(
            splicing_event_count.items(), key=lambda x: x[0]
        )
        for event, count in splicing_count_sorted:
            chrom, start, end = event
            # convert to 1-based coordinates
            out_f.write(f"{chrom}\t{start + 1}\t{end}\t{count}\n")
