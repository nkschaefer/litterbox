#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from cleanup_cat_annotation import get_tags, join_tags
import subprocess
import random
from collections import Counter
"""
Given a list of genes that were removed from the CAT annotation because all transcripts
were filtered out, a mapping of these (human) genes to Ensembl gene IDs for the species of
interest, an old Ensembl gene annotation for the species of interest on the same genome 
assembly, attempts to rescue these genes by taking their records from the old Ensembl 
annotation.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf", 
        "-g", 
        help="Ensembl annotation for the target (non-human) species, \
on the same assembly as the CAT annotation. Can be gzipped.",
        required=True
    )
    parser.add_argument(
        "--genelist", 
        "-G", 
        help="List of genes from get_homologous_genelist.R. These \
should be human<tab>otherspecies Ensembl IDs for genes that were \
discarded from the CAT annotation for having all unsuitable \
transcripts.", 
        required=True
    )
    parser.add_argument(
        "--output_prefix", 
        "-o", 
        required=True
    )
    parser.add_argument(
        "--hgnc_ens", 
        "-he", 
        help="File with 3 (tab separated) columns: HGNC ID, HGNC \
approved symbol, Ensembl ID (or blank if missing)", 
        required=True
    )
    parser.add_argument(
        "--idmap",
        help="File mapping from Ensembl seq ID to CAT seq ID",
        required=True
    )
    return parser.parse_args()

def main(args):
    options = parse_args()
    
    ens2name = {}
    hgnc2name = {}

    # Read HGNC file
    f = open(options.hgnc_ens, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if len(dat) > 2 and len(dat[2]) > 0 and len(dat[1]) > 0:
            ens2name[dat[2]] = dat[1]
        if len(dat[1]) > 0:
            hgnc2name[dat[0]] = dat[1]
    f.close()
    
    gene2human = {}
    human2gene = {}
    humancounts = Counter
    f = open(options.genelist, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        gene2human[dat[1]] = dat[0]
        human2gene[dat[0]] = dat[1]
        humancounts[dat[0]] += 1
    f.close()
    # Remove human genes that map to more than 1 other species gene
    for human in human2gene:
        other = human2gene[human]
        if humancounts[human] > 1:
            del gene2human[other]

    tx2gene = {}

    outf = open("{}.rescued.gtf".format(options.output_prefix), 'w')
    out_list = open("{}.rescued.gid".format(options.output_prefix), 'w')
    
    id2id = {}
    f = open(options.idmap, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            dat = line.split('\t')
            id2id[dat[0]] = dat[1]
    f.close()

    is_gz = False
    if options.gtf == '-':
        f = sys.stdin
    elif options.gtf[-3:] == '.gz':
        is_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f = open(options.gtf, 'r')
    for line in f:
        
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        
        if line[0] == '#':
            continue
        dat = line.split('\t')
        
        dat[0] = id2id[dat[0]]
        
        # Shouldn't happen but check anyway
        if dat[0] == 'chrM' or dat[0] == 'MT' or dat[0] == 'chrMT':
            continue
        
        tags = get_tags(dat)
        if 'gene_id' in tags:
            if tags['gene_id'] in gene2human:
                ensg = gene2human[tags['gene_id']]
                if ensg in ens2name:
                    tags['gene_name'] = ens2name[ensg]
                else:
                    tags['gene_name'] = ensg
                tags['source_gene'] = ensg
                print(ensg, file=out_list)

                # No need to alter gene_id or transcript_id
                # Tell users which annotation we took this from
                dat[2] = options.gtf.split('/')[-1]

                if 'transcript_id' in tags:
                    tx2gene[tags['transcript_id']] = ensg
                
                dat[8] = join_tags(tags)
                print("\t".join(dat), file=outf)
    
    if options.gtf != '-':
        f.close()
    outf.close()
    out_list.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
