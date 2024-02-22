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
interest, an old Ensembl gene annotation for the species of interest, and a UCSC chain 
file for mapping from the old genome assembly to the new, attempts to rescue these genes 
by lifting over their records from the old Ensembl annotation to the new genome assembly.
Requires gtfToGenePred, liftOver, and genePredToGtf from the UCSC Kent utilities.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf", 
        "-g", 
        help="Ensembl annotation on an older assembly of the target \
(non-human) species. Can be gzipped.", 
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
        "--chain", 
        "-c", 
        help="Chain file for liftOver", 
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
        help="File mapping UCSC sequence ID to sequence IDs in the \
HAL-aligned assembly version on which CAT was run",
        required=True
    )
    parser.add_argument(
        "--gtfToGenePred_path", 
        help="Path to gtfToGenePred (OPTIONAL; default = in path)",
        required=False
    )
    parser.add_argument(
        "--liftOver_path", 
        help="Path to liftOver (OPTIONAL; default = in path)", 
        required=False
    )
    parser.add_argument(
        "--genePredToGtf_path",
        help="Path to genePredToGtf (OPTIONAL; default = in path)",
        required=False
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
    human_counts = Counter()
    f = open(options.genelist, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        gene2human[dat[1]] = dat[0]
        human2gene[dat[0]] = dat[1]
        human_counts[dat[0]] += 1
    f.close()

    # Toss out human genes that map to more than 1 gene in other
    # species
    for human in human2gene:
        other = human2gene[human]
        if human_counts[human] > 1:
            del gene2human[other]
    
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
        
        if is_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        
        if line[0] == '#':
            continue
        dat = line.split('\t')
        
        # Here, we are assuming that we're using the UCSC liftOver files to go from
        # one UCSC genome assembly to another, but we're using an Ensembl annotation.
        # That means we're likely to need to append "chr" to the beginning of 
        # chromosome names.
        if len(dat[0]) < 4 or dat[0:3] != "chr":
            dat[0] = "chr" + dat[0]
        if dat[0] == "chrM" or dat[0] == "chrMT":
            continue
        
        tags = get_tags(dat)
        if 'gene_id' in tags:
            if tags['gene_id'] in gene2human:
                ensg = gene2human[tags['gene_id']]
                tags['gene_id'] = gene2human[tags['gene_id']]
                if ensg in ens2name:
                    tags['gene_name'] = ens2name[ensg]
                else:
                    tags['gene_name'] = ensg
                if 'transcript_id' in tags:
                    tx2gene[tags['transcript_id']] = ensg
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

    subprocess.call([
        gtfToGenePred, 
        '-genePredExt', 
        '{}.old.gtf'.format(options.output_prefix),
        '{}.old.gp'.format(options.output_prefix)
    ])
    subprocess.call([
        liftOver, 
        '-genePred', 
        '{}.old.gp'.format(options.output_prefix), 
        options.chain, 
        '{}.new.gp'.format(options.output_prefix), 
        '{}.rescue.unmapped'.format(options.output_prefix)
    ])
    subprocess.call([
        genePredToGtf, 
        '-utr', 
        'file', 
        '{}.new.gp'.format(options.output_prefix),
        '{}.new.gtf'.format(options.output_prefix)
    ])
    
    # Load ID mappings
    id2id = {}
    f = open(options.idmap, 'r')
    for line in f:
        line = line.rstrip()
        if line != "":
            dat = line.split('\t')
            id2id[dat[0]] = dat[1]
    f.close()

    f = open('{}.new.gtf'.format(options.output_prefix), 'r')
    f_out = open('{}.rescued.gtf'.format(options.output_prefix), 'w')
    f_out_list = open('{}.rescued.gid'.format(options.output_prefix), 'w')

    rand_tag = ''.join(random.choice('0123456789ABCDEF') for i in range(5))
    
    # Make two passes: first pass is to learn what gene IDs have a gene feature
    # that made it through liftOver. Second is to actually print stuff
    has_gene = set([])
    has_tx = set([])
    has_exon = set([])
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        tags = get_tags(dat)
        if dat[2] == 'gene':
            has_gene.add(tags['gene_id'])
        elif dat[2] == 'transcript':
            has_tx.add(tags['gene_id'])
        elif dat[2] == 'exon':
            has_exon.add(tags['gene_id'])
    f.close()
    
    f = open('{}.new.gtf'.format(options.output_prefix), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        # Convert UCSC seq ID to HAL-aligned assembly sequence ID
        dat[0] = id2id[dat[0]]
        # Mark source
        dat[1] = 'liftOver_{}'.format(options.gtf.split('/')[-1])
        tags = get_tags(dat)
        # Make sure this gene has a gene entry plus at least one transcript
        # and at least one exon.
        if tags['gene_id'] in has_gene and tags['gene_id'] in has_tx and \
            tags['gene_id'] in has_exon:
            if tags['gene_id'] in ens2name:
                tags['gene_name'] = ens2name[tags['gene_id']]
            tags['source_gene'] = tags['gene_id']
            print(tags['gene_id'], file=f_out_list)
            # Make gene IDs unique. Probably not important, but could potentially have a collision - 
            # we want to note that this gene is not the same thing as the source human gene
            tags['gene_id'] += '-' + rand_tag
            dat[8] = join_tags(tags)
            print("\t".join(dat), file=f_out)
    f_out.close()
    f_out_list.close()
    f.close()
    
    # Clean up
    os.unlink('{}.old.gtf'.format(options.output_prefix))
    os.unlink('{}.old.gp'.format(options.output_prefix))
    os.unlink('{}.new.gp'.format(options.output_prefix))
    os.unlink('{}.new.gtf'.format(options.output_prefix))


if __name__ == '__main__':
    sys.exit(main(sys.argv))
