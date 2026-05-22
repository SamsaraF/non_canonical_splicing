#!/usr/bin/env python3
# -*- coding=utf-8 -*-

from argparse import ArgumentParser


def main(bam_list, output_plot, coord=None, bam_labels=None):
    import subprocess
    import seaborn as sns
    import pandas as pd
    import io

    cmd = f'samtools depth -a {" ".join(bam_list)}'
    if coord:
        cmd += f' -r {coord}'
    try:
        depth = subprocess.run(cmd, capture_output=True, text=True, shell=True).stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools depth: {e}")
        return
    
    names = ['chrom', 'pos']
    if bam_labels:
        names += bam_labels
    else:
        names += [i for i in range(len(bam_list))]
    depth_df = pd.read_table(io.StringIO(depth), header=None, names=names)
    plot_df = depth_df.melt(id_vars=['chrom', 'pos'], var_name='sample', value_name='depth')

    sns.set_style('whitegrid')
    for chrom, sub_df in plot_df.groupby('chrom'):
        plot = sns.lineplot(data=sub_df, x='pos', y='depth', hue='sample')
        plot.set_title(chrom)
        plot.set_xlabel('Position')
        plot.set_ylabel('Depth')
        plot.get_figure().savefig(f"{chrom}.{output_plot}", dpi=300)
        plot.get_figure().clf()


if __name__ == '__main__':
    parser = ArgumentParser(description='Visualize samtools depth output')
    parser.add_argument('-b', '--bam', nargs='+', required=True, help='Path to the BAM files. Can accept multiple files for comparison.')
    parser.add_argument('-l', '--label', nargs='+', help='Labels for the BAM files, in the same order as the BAM files. If not provided, index (0, 1, 2, ...) will be used as labels.')
    parser.add_argument('-o', '--output', required=True, help='Path to save the output plot')
    parser.add_argument('-c', '--coord', help='coordinate to focus on, in the format chrom:start-end (e.g., chr1:1000-2000). If not provided, the entire depth will be plotted.')
    args = parser.parse_args()

    if args.label and len(args.label) != len(args.bam):
        print("Error: The number of labels must match the number of BAM files.")
        exit(1)

    main(args.bam, args.output, args.coord, args.label)
