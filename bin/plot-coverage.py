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


def parse_depths(pileup_path: Path) -> list[dict]:
    """
    Parse pileup file

    :param pileup_path: Pileup file
    :type pileup_path: Path
    :return: List of dictionaries with the keys chrom, pos, ref, and depth
    :rtype: list[dict]
    """
    depths = []
    num_positions = 0
    previous_chrom = None
    with open(pileup_path, 'r') as f:
        for line in f:
            num_positions += 1
            line_split = line.strip().split("\t")
            chrom = line_split[0]
            pos = int(line_split[1])
            depth = int(line_split[2])
            output = {
                "chrom": chrom,
                "pos": pos,
                "depth": depth
            }
            depths.append(output)
            if previous_chrom is None:
                previous_chrom = chrom
            elif previous_chrom != chrom:
                num_positions = 0
                previous_chrom = chrom
            if num_positions % 100_000 == 0:
                logging.info(f"Processed {num_positions:,} positions from {chrom}...")

    return depths


def plot_coverage(depths: pd.DataFrame, bed: dict, sample_name: str="Sample", threshold: int=10, y_limit: int=500, log_scale: bool=False, width_inches_per_mb=3, height_inches_per_chrom=2) -> sns.FacetGrid:
    """
    Plot coverage

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
    :param width_inches_per_mb: Width inches per megabase (default: 12)
    :type width_inches_per_mb: float
    :param height_inches_per_chrom: Height inches per chromosome (default: 8)
    :type height_inches_per_chrom: float
    :return: Seaborn FacetGrid object
    :rtype: seaborn.FacetGrid
    """
    try:
        percent_coverage_above_threshold = (
            sum(1 if x > threshold else 0 for x in depths.depth)
            / depths.shape[0]
            * 100
        )
        logging.info(f"Percent bases with coverage above {threshold}X: {percent_coverage_above_threshold: .1f}%")
    except ZeroDivisionError:
        percent_coverage_above_threshold = 0

    # Filter out any chroms with fewer than 100 positions
    chrom_counts = depths.chrom.value_counts()
    chroms_to_keep = chrom_counts[chrom_counts >= 100].index
    depths = depths[depths.chrom.isin(chroms_to_keep)]
    num_chroms = depths.chrom.nunique()
    longest_chrom = depths.chrom.value_counts().idxmax()
    longest_chrom_length = depths[depths.chrom == longest_chrom].pos.max()
    longest_chrom_length_mb = longest_chrom_length / 1_000_000
    sns.set_style("whitegrid")
    dpi = 300
    width_pixels = round(longest_chrom_length_mb * width_inches_per_mb * dpi)
    height_pixels = round(height_inches_per_chrom * num_chroms * dpi)
    logging.info(f"Plot size: {width_pixels} x {height_pixels} pixels")
    try:
        g = sns.FacetGrid(depths, row="chrom", height=height_inches_per_chrom, aspect=(width_inches_per_mb * longest_chrom_length_mb)/height_inches_per_chrom, gridspec_kws={"hspace": 0.5})
    except ValueError as e:
        logging.error(e)
        logging.error("Try reducing the width_inches_per_mb and height_inches_per_chrom parameters.")
        exit(1)
    g.fig.set_size_inches(longest_chrom_length_mb * width_inches_per_mb, height_inches_per_chrom * num_chroms)
    g.map(sns.lineplot, "pos", "depth", linewidth=0.5)
    g.map(plt.axhline, y=threshold, color="red", linestyle="--", linewidth=0.5)
    g.set_xlabels("Position")
    for ax in g.axes.flatten():
        ax.tick_params(labelbottom=True)
    g.set(xlim=(0, longest_chrom_length))
    if log_scale:
        g.set(yscale="log")
        g.set(ylim=(1, y_limit))
        g.set_ylabels("Depth of coverage (log scale)")
    else:
        g.set(ylim=(0, y_limit))
        g.set_ylabels("Depth of coverage")
    # set titles for each facet with chrom name
    g.set_titles(row_template="{row_name}")

    # line_plot.set_xticks(range(0, depths.shape[0], 100000), minor=True)
    
    plt.close()

    return g



def main(args):

    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        level=logging.INFO,
        datefmt='%Y-%m-%dT%H:%M:%S'
    )

    all_depths_df = pd.read_csv(args.depths, sep="\t")
    # Check if depths are empty
    if all_depths_df.empty:
        logging.error("No depth of coverage data found in input file. Exiting...")
        exit(0)

    logging.info("Calculating tumbling window depths")
    # create an empty dataframe to fill with tumbling window depths:
    tumbling_window_depths = pd.DataFrame(columns=["chrom", "pos", "depth"])
    for chrom, group in all_depths_df.groupby("chrom"):
        # for each chromosome, iterate over the depths in tumbling windows:
        for i in range(0, len(group), args.window):
            # get a slice of the dataframe for the tumbling window:
            window = group.iloc[i:i+args.window]
            # calculate the middle position of the window:
            middle_pos = window.iloc[int(len(window)/2)].pos
            # calculate the mean depth of the window:
            mean_depth = window.depth.mean()
            record = {
                "chrom": chrom,
                "pos": middle_pos,
                "depth": mean_depth
            }
            # Avoid FutureWarning about empty or all-NA entries:
            if len(tumbling_window_depths) == 0:
                tumbling_window_depths = pd.DataFrame([record])
            else:
                tumbling_window_depths = pd.concat([tumbling_window_depths, pd.DataFrame([record])], ignore_index=True)
    tumbling_window_depths = tumbling_window_depths.reset_index(drop=True)

    bed = {}
    if args.bed is not None:
        bed = parse_bed(args.bed)

    logging.info("Plotting coverage...")

    coverage_plot = plot_coverage(tumbling_window_depths, bed, args.sample_id, args.threshold, args.y_limit, args.log_scale, args.width_inches_per_mb, args.height_inches_per_chrom)
    logging.info(f"Saving plot to {args.output}...")
    coverage_plot.savefig(args.output, dpi=100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot coverage')
    parser.add_argument('-d', '--depths', help='Input .bam file')
    parser.add_argument('-r', '--ref', help='Input ref .fasta file', required=True)
    parser.add_argument('--bed', type=Path, help='bed file (optional)')
    parser.add_argument('-s', '--sample-id', default="Sample", help='Sample ID')
    parser.add_argument('-t', '--threshold', type=int, default=10, help='Threshold for depth of coverage')
    parser.add_argument('--window', type=int, default=1000, help='Tumbling window for coverage')
    parser.add_argument('--log-scale', action='store_true', help='Log scale y-axis')
    parser.add_argument('--width-inches-per-mb', type=float, default=3, help='Width inches per megabase (default: 3)')
    parser.add_argument('--height-inches-per-chrom', type=float, default=2, help='Height inches per chromosome (default: 2)')
    parser.add_argument('--y-limit', type=int, default=500, help='Maximum y-axis value in plot (default: 500)')
    parser.add_argument('-o', '--output', default="./coverage.png", help='Output file')
    args = parser.parse_args()
    main(args)
