#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from cleanup_cat_annotation import get_tags, join_tags
import subprocess
"""
Given a list of genes that were removed from the CAT annotation because all transcripts
were filtered out, a mapping of these (human) genes to Ensembl gene IDs for the species of
interest, an old Ensembl gene annotation for the species of interest, and a UCSC chain 
file for mapping from the old genome assembly to the new, attempts to rescue these genes
by lifting over their records from the old Ensembl annotation to the new genome assembly.
Requires gtfToGenePred, liftOver, and genePredToGtf from the UCSC Kent utilities to be
present in the same directory where this is run.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", "-g", help="Ensembl annotation on an older assembly of the \
target (non-human) species. Can be gzipped.", required=True)
    parser.add_argument("--genelist", "-G", help="List of genes from get_homologous_genelist.R. \
These should be human<tab>otherspecies Ensembl IDs for genes that were discarded from the \
CAT annotation for having all unsuitable transcripts.", required=True)
    parser.add_argument("--output_prefix", "-o", required=True)
    parser.add_argument("--chain", "-c", help="Chain file for liftOver", required=True)
    parser.add_argument("--hgnc_ens", "-he", help="File with 3 (tab separated) columns: \
HGNC ID, HGNC approved symbol, Ensembl ID (or blank if missing)", required=True)
    parser.add_argument("--gtfToGenePred_path", help="Path to gtfToGenePred (OPTIONAL; default = in path)", \
        required=False)
    parser.add_argument("--liftOver_path", help="Path to liftOver (OPTIONAL; default = in path)", required=False)
    parser.add_argument("--genePredToGtf_path", help="Path to genePredToGtf (OPTIONAL; default = in path)", \
        required=False)
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
    f = open(options.genelist, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        gene2human[dat[1]] = dat[0]
    f.close()
    
    tx2gene = {}

    outf = open("{}.old.gtf".format(options.output_prefix), 'w')

    is_gz = False
    if options.gtf == '-':
        f = sys.stdin
    elif options.gtf[-3:] == '.gz':
        is_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f = open(options.gtf, 'r')
    for line in f:
        if line[0] == '#':
            continue
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        dat = line.split('\t')
        tags = get_tags(dat)
        if 'gene_id' in tags:
            if tags['gene_id'] in gene2human:
                tags['gene_id'] = gene2human[tags['gene_id']]
                if tags['gene_id'] in ens2name:
                    tags['gene_name'] = ens2name[tags['gene_id']]
                if 'transcript_id' in tags:
                    tx2gene[tags['transcript_id']] = tags['gene_id']
                dat[8] = join_tags(tags)
                print("\t".join(dat), file=outf)
    
    if options.gtf != '-':
        f.close()
    outf.close()

    gtfToGenePred = 'gtfToGenePred'
    if options.gtfToGenePred_path is not None:
        if options.gtfToGenePred_path[-1] == '/':
            gtfToGenePred = options.gtfToGenePred_path + 'gtfToGenePred'
        else:
            gtfToGenePred = options.gtfToGenePred_path + '/gtfToGenePred'
    liftOver = 'liftOver'
    if options.liftOver_path is not None:
        if options.liftOver_path[-1] == '/':
            liftOver = options.liftOver_path + 'liftOver'
        else:
            liftOver = options.liftOver_path + '/liftOver'
    genePredToGtf = 'genePredToGtf'
    if options.genePredToGtf_path is not None:
        if options.genePredToGtf_path[-1] == '/':
            genePredToGtf = options.genePredToGtf_path + 'genePredToGtf'
        else:
            genePredToGtf = options.genePredToGtf_path + '/genePredToGtf'

    subprocess.call([gtfToGenePred, '-genePredExt', '{}.old.gtf'.format(options.output_prefix), \
        '{}.old.gp'.format(options.output_prefix)])
    subprocess.call([liftOver, '-genePred', '{}.old.gp'.format(options.output_prefix), \
        options.chain, '{}.new.gp'.format(options.output_prefix), '{}.unmapped'.format(options.output_prefix)])
    subprocess.call([genePredToGtf, '-utr', 'file', '{}.new.gp'.format(options.output_prefix), \
        '{}.new.gtf'.format(options.output_prefix)])
    
    f = open('{}.new.gtf'.format(options.output_prefix), 'r')
    f_out = open('{}.rescued.gtf'.format(options.output_prefix), 'w')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        tags = get_tags(dat)
        if tags['gene_id'] in ens2name:
            tags['gene_name'] = ens2name[tags['gene_id']]
        dat[8] = join_tags(tags)
        print("\t".join(dat), file=f_out)
    f_out.close()
    f.close()
    
    # Clean up
    os.unlink('{}.old.gtf'.format(options.output_prefix))
    os.unlink('{}.old.gp'.format(options.output_prefix))
    os.unlink('{}.new.gp'.format(options.output_prefix))
    os.unlink('{}.new.gtf'.format(options.output_prefix))


if __name__ == '__main__':
    sys.exit(main(sys.argv))
