#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from collections import Counter, defaultdict
from cleanup_cat_annotation import get_tags, join_tags 
"""
After filtering the CAT annotation, put the gene symbols
chosen to represent the genes (i.e. from HUGO) into the
human annotation (i.e. from 10X).

Also discard any genes that were filtered out.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf", 
        "-g", 
        help="Input GTF file for human. Can be gzipped.", 
        required=True
    )
    parser.add_argument(
        "--allowed_genes", 
        "-G", 
        help="List of allowed gene IDs and their names.", 
        required=True
    )
    parser.add_argument(
        "--allowed_tx",
        "-T",
        help="List of allowed Ensembl transcript IDs",
        required=True
    )
    parser.add_argument(
        "--tx2gene",
        "-t2g",
        help="Map of transcript ID to gene ID (tab separated, all tx/gene combos)",
        required=True
    )
    parser.add_argument(
        "--mito", 
        "-m", 
        help="List of allowed gene IDs and and their names, \
on the mitochondrial genome.", 
        required=False
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output GTF file to create",
        required=True
    )
    parser.add_argument(
        "--mito_name",
        "-M",
        help="Name of mitochondrial genome in the human annotation (OPTIONAL, default chrM)",
        required=False
    )
    return parser.parse_args()

"""
Check whether a file is gzipped.
"""
def file_is_gz(file):
    with open(file, 'rb') as test:
        return test.read(2) == b'\x1f\x8b'


def main(args):
    options = parse_args()
    
    # First, load which transcripts are from which genes
    tx2gene = {}
    f = open(options.tx2gene, 'r')
    for line in f:
        line = line.rstrip()
        if len(line) > 0:
            dat = line.split('\t')
            tx2gene[dat[0]] = dat[1]
    f.close()

    # Next, load all allowed transcripts and store counts of allowed
    # transcripts per gene
    tx_allowed = set([])
    gene_tx_counts = Counter()
    f = open(options.allowed_tx, 'r')
    for line in f:
        line = line.rstrip()
        if len(line) > 0:
            tx_allowed.add(line)
            gene_tx_counts[tx2gene[line]] += 1
    f.close()
    
    # Finally, load allowed genes and remove any that have 0 valid
    # transcripts
    ensg2n = {}
    f = open(options.allowed_genes, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if len(dat) >= 2:
            ensg = dat[0]
            name = dat[1]
            if ensg in gene_tx_counts and gene_tx_counts[ensg] > 0:
                ensg2n[ensg] = name
    f.close()
    
    if options.mito is not None:
        f = open(options.mito, 'r')
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            if len(dat) >= 2:
                # Assume all mitochondrial genes are valid
                ensg2n[dat[0]] = dat[1]
        f.close()

    outf = open(options.output, 'w')

    f = None
    is_gz = file_is_gz(options.gtf)
    if is_gz:
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
            par = False 
            if dat[0] == options.mito_name:
                # Ensure we keep all records on the mitochondrion
                print(line, file=outf)
            else:
                if 'gene_id' in tags:
                    if '_PAR_Y' in tags['gene_id']:
                        # Skip this
                        par = True
                    else:
                        ensg = tags['gene_id'].split('.')[0]
                        if ensg in ensg2n:
                            tags['gene_name'] = ensg2n[ensg]
                            tags['gene_id'] = ensg
                            print_feature = True
                            if 'transcript_id' in tags:
                                # Also ensure it's a valid transcript.
                                txid = tags['transcript_id'].split('.')[0]
                                if txid in tx_allowed:
                                    print_feature = True
                                else:
                                    print_feature = False
                                tags['transcript_id'] = txid
                            if 'tag' in tags:
                                for x in tags['tag']:
                                    if x == 'readthrough_transcript' or x == 'PAR' or \
                                        x == 'stop_codon_readthrough' or x == 'low_sequence_quality':
                                        print_feature = False
                                        break
                            if print_feature:
                                dat[8] = join_tags(tags)
                                print("\t".join(dat), file=outf)
    outf.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

