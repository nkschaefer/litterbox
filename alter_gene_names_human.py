#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from cleanup_cat_annotation import get_tags, join_tags 
"""
After filtering the CAT annotation, put the gene symbols
chosen to represent the genes (i.e. from HUGO) into the
human annotation (i.e. from 10X).
"""

def parse_args():
    parser = argparse.ArgumentParser(description='__doc__')
    parser.add_argument("--gtf", "-g", help="Input GTF file for human. \
Can be gzipped.", required=True)
    parser.add_argument("--allowlist", "-a", help="List of allowed gene IDs \
and their names.", required=True)
    parser.add_argument("--mito", "-m", help="List of allowed gene IDs and \
and their names, on the mitochondrial genome.", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    ensg2n = {}
    f = open(options.allowlist, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if len(dat) >= 2:
            ensg2n[dat[0]] = dat[1]
    f.close()
    if options.mito is not None:
        f = open(options.mito, 'r')
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            if len(dat) >= 2:
                ensg2n[dat[0]] = dat[1]
        f.close()

    f = None
    is_gz = False
    if options.gtf[-3:] == '.gz':
        is_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f = open(options.gtf, 'r')

    for line in f:
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#" and line != "":
            dat = line.split('\t')
            tags = get_tags(dat)
            if 'gene_id' in tags and tags['gene_id'] in ensg2n:
                tags['gene_name'] = ensg2n[tags['gene_id']]
                dat[8] = join_tags(tags)
            print("\t".join(dat))


if __name__ == '__main__':
    sys.exit(main(sys.argv))

