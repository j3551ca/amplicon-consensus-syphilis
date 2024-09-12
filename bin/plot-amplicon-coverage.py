#!/usr/bin/env python3

import argparse
import csv
import json
import logging
import os
import subprocess

from pathlib import Path
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import pandas as pd
import seaborn as sns

from matplotlib.patches import Rectangle

def parse_bed(bed_path: Path) -> dict:
    """
    Parse BED file

    :param bed: BED file
    :type bed: Path
    :return: Dictionary of genes by name
    :rtype: dict
    """
    genes_by_name = {}
    with open(bed_path, 'r') as f:
        for line in f:
            line_split = line.strip().split('\t')
            chrom = line_split[0]
            start = int(line_split[1])
            end = int(line_split[2])
            name = line_split[3]
            strand = line_split[5]
            genes_by_name[name] = {
                'name': name,
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand,
            }

    return genes_by_name


def parse_amplicon_depths(amplicon_depths_path: Path) -> list[dict]:
    """
    Parse a tsv file with amplicon depths

    :param amplicon_depths_path: Pileup file
    :type amplicon_depths_path: Path
    :return: List of dictionaries with the keys
    :rtype: list[dict]
    """
    depths = []
    
    with open(amplicon_depths_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            record = {
                'reference_name': row['reference_name'],
                'amplicon_id': int(row['amplicon_id']),
                'mean_depth': int(row['mean_depth']),
            }
            depths.append(record)

    return depths


def plot_coverage(depths: pd.DataFrame, sample_name: str="Sample", threshold: int=10, y_limit: int=500, log_scale: bool=False) -> matplotlib.figure.Figure:
    """
    Create bar plot of coverage for each amplicon

    :param depths: DataFrame of depths (columns: chrom, pos, depth)
    :type depths: pandas.DataFrame
    :param bed: Dictionary of bed regions by name
    :type bed: dict
    :param sample_name: Sample name
    :type sample_name: str
    :param threshold: Threshold for depth of coverage
    :type threshold: int
    :param y_limit: Maximum y-axis value in plot (default: 500)
    :type y_limit: int
    :param log_scale: Log scale y-axis (default: False)
    :type log_scale: bool
    :return: The figure
    :rtype: matplotlib.figure.Figure
    """
    sns.set_theme()
    sns.set_style("whitegrid")
    
    ax = sns.barplot(
        data=depths,
        x=(depths.index + 1),
        y='mean_depth',
    )
        

    ax.set(
        xlabel='Amplicon',
        ylabel='Mean Depth',
        title=f'{sample_name} Amplicon Coverage',
    )

    if log_scale:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(mtick.ScalarFormatter())
        ax.yaxis.set_minor_formatter(mtick.NullFormatter())
    else:
        ax.set_ylim(0, y_limit)

    # Add a horizontal line at the threshold
    ax.axhline(y=threshold, color='r', linestyle='--')
    
    fig = ax.get_figure()

    # Set the size of the plot based on the number of amplicons
    num_amplicons = depths.shape[0]
    fig.set_size_inches(0.5 * num_amplicons, 6)
    
    return fig



def main(args):

    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        level=logging.INFO,
        datefmt='%Y-%m-%dT%H:%M:%S'
    )

    amplicon_depths_df = pd.read_csv(args.depths, sep="\t")

    if amplicon_depths_df.empty:
        logging.error("No depth of coverage data found in input file. Exiting...")
        exit(0)

    coverage_plot = plot_coverage(amplicon_depths_df, sample_name=args.sample_id, threshold=args.threshold, y_limit=args.y_limit, log_scale=args.log_scale)
    coverage_plot.savefig(args.output, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coverage')
    parser.add_argument('-d', '--depths', help='Input .bam file')
    parser.add_argument('-s', '--sample-id', default="Sample", help='Sample ID')
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    parser.add_argument('--log-scale', action='store_true', help='Log scale y-axis')
    parser.add_argument('--y-limit', type=int, default=2000, help='Maximum y-axis value in plot (default: 2000)')
    parser.add_argument('-o', '--output', default="./coverage.png", help='Output file')
    args = parser.parse_args()
    main(args)
