#!/usr/bin/env python3
import sys
import os
import argparse
import csv
import natsort
import pandas
"""
This code comes from here:
    https://github.com/10XGenomics/cellranger/issues/133
And was written by chbk:
    https://github.com/chbk

I adapted it slightly to take in (and concatenate) multiple input GTFs.

It is desigend to solve an issue where cellranger-arc mkref complains about
the sort order of a GTF. Grouped records need to be near each other in the
file and position sorting will not work.
"""

def parse_args():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", "-i",
        required = True,
        help = 'Input GTFs',
        nargs="+",
    )
    parser.add_argument(
        "--output", "-o",
        required = True,
        help = 'Output sorted gtf'
    )
    args = parser.parse_args()
    return args

def main(args):
    options = parse_args()

    # Read gtf input
    gtf_columns = {
        'chromosome': 'str',
        'source': 'str',
        'feature': 'str',
        'start': 'uint64',
        'end': 'uint64',
        'score': 'str',
        'strand': 'str',
        'frame': 'str',
        'attribute': 'str'
    }

    gtfs = []
    for name in options.input:
        gtfs.append(pandas.read_csv(
          name,
          sep = '\t',
          comment = '#',
          names = gtf_columns.keys(),
          dtype = gtf_columns
        ))
    
    gtf = pandas.concat(gtfs, axis=0, ignore_index=True)
    
    # Get gene id and transcript id for each row
    gtf['gene'] = gtf['attribute'].str.extract(r'gene_id "(.+?)"')
    gtf['transcript'] = gtf['attribute'].str.extract(r'transcript_id "(.+?)"')

    # Get genes start and end positions
    genes = gtf[gtf['feature'] == 'gene'][
        ['chromosome', 'strand', 'gene', 'start', 'end']
    ].rename(
        columns = {
            'start': 'gene_start',
            'end': 'gene_end'
        }
    ).set_index(
        ['chromosome', 'strand', 'gene']
    )

    # Get transcripts start and end positions
    transcripts = gtf[gtf['feature'] == 'transcript'][
        ['chromosome', 'strand', 'transcript', 'start', 'end']
    ].rename(
        columns = {
            'start': 'transcript_start',
            'end': 'transcript_end',
        }
    ).set_index(
        ['chromosome', 'strand', 'transcript']
    )

    # Add gene and transcript start and end positions to each row
    gtf = gtf.set_index(
        ['chromosome', 'strand', 'gene']
    ).merge(
        genes,
        how = 'left',
        on = ['chromosome', 'strand', 'gene']
    ).reset_index().set_index(
        ['chromosome', 'strand', 'transcript']
    ).merge(
        transcripts,
        how = 'left',
        on = ['chromosome', 'strand', 'transcript']
    ).reset_index()

    # Sort rows
    gtf = gtf.sort_values(
        by = [
            'chromosome',
            'strand',
            'gene_start',
            'gene_end',
            'gene',
            'transcript_start',
            'transcript_end',
            'transcript',
            'feature',
            'start',
            'end'
        ],
        key = lambda x: (
            [0 if i == 'gene' else 1 if i == 'transcript' else 2 for i in x]
            if x.name == 'feature'
            else natsort.natsort_key(x)
        )
    )

    # Write gtf to output
    gtf.to_csv(
        options.output,
        sep = '\t',
        columns = gtf_columns.keys(),
        header = False,
        index = False,
        quoting = csv.QUOTE_NONE,
        float_format = '%.10g'
    )

if __name__ == '__main__':
    sys.exit(main(sys.argv))

