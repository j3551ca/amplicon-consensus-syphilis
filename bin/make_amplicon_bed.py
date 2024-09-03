#!/usr/bin/env python

import argparse
import csv
import json
import os
import sys


def parse_primer_bed(primer_bed):
    """
    Parse a BED file with primer coordinates.

    :param primer_bed: Path to the BED file with primer coordinates
    :type primer_bed: str
    :return: List of dictionaries with primer coordinates
    :rtype: list
    """
    primers = []
    # check if the file has a header
    has_header = False
    with open(primer_bed, 'r') as f:
        header = f.readline()
        if header.startswith('chrom'):
            has_header = True

    if has_header:
        with open(primer_bed, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                primer = {
                    'chrom': row['chrom'],
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'name': row['name'],
                }
                primers.append(primer)
    else:
        with open(primer_bed, 'r') as f:
            for line in f:
                line_split = line.strip().split('\t')
                primer = {
                    'chrom': line_split[0],
                    'start': int(line_split[1]),
                    'end': int(line_split[2]),
                    'name': line_split[3],
                }
                primers.append(primer)

    return primers


def make_amplicon_bed_record(primer1, primer2):
    """
    Make an amplicon BED record from two primers.

    The amplicon is the region between the end of primer1 and the start of primer2.

    We check that the primers are on the same chromosome and that they do not overlap.

    :param primer1: Dictionary with primer coordinates
    :type primer1: dict
    :param primer2: Dictionary with primer coordinates
    :type primer2: dict
    :return: Dictionary with amplicon coordinates
    :rtype: dict
    """
    chrom = primer1['chrom']
    primer1_end = primer1['end']
    primer2_start = primer2['start']

    if primer1['chrom'] != primer2['chrom']:
        sys.exit(f"Error: The primers are on different chromosomes: {primer1['chrom']} and {primer2['chrom']}")

    if primer1_end > primer2_start:
        sys.exit(f"Error: The primers overlap: {primer1_end} and {primer2_start}")

    amplicon = {
        'chrom': chrom,
        'start': primer1_end,
        'end': primer2_start,
        'name': f"{primer1['name']}_{primer2['name']}",
    }

    return amplicon


def main(args):
    primers = parse_primer_bed(args.primer_bed)

    # check there are an even number of primers
    if len(primers) % 2 != 0:
        sys.exit(f"Error: The number of primers in the BED file is not even: {len(primers)}")

    # We assume that the primers are in pairs, with the first primer in the pair being the forward primer
    # and the second primer in the pair being the reverse primer.
    output = []
    for i in range(0, len(primers), 2):
        primer1 = primers[i]
        primer2 = primers[i + 1]
        amplicon = make_amplicon_bed_record(primer1, primer2)
        output.append(amplicon)

    output_fieldnames = [
        'chrom',
        'start',
        'end',
        'name'
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='excel-tab', quoting=csv.QUOTE_MINIMAL, extrasaction='ignore')
    writer.writeheader()

    for record in output:
        writer.writerow(record)
        
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Description of your script')
    parser.add_argument('--primer-bed', type=str, required=True, help='BED file with primer coordinates')
    args = parser.parse_args()
    main(args)
