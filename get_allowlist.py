#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
from cleanup_cat_annotation import join_tags
"""
Given the most recent GENCODE annotation GTF, finds all genes that would
pass filter (according to 10X Genomics guidelines).

Also reads in preferred gene common names from HGNC. Prints a list of 
genes passing filter as a tab separated file of Ensembl gene ID, accepted name.

Finally, can print a list of all genes on the mitochondrial genome, so they
can be excluded from the CAT annotation.
"""
def parse_args():
    parser = argparse.ArgumentParser(description='__doc__')
    parser.add_argument("--gtf", "-g", \
        help="Latest GENCODE gtf (can be gzipped)", required=True)
    parser.add_argument("--hgnc_ens", "-he", help="File with 3 (tab separated) columns: \
HGNC ID, HGNC approved symbol, Ensembl ID (or blank if missing)", required=True)
    parser.add_argument("--output_prefix", "-o", help="Prefix for output files. Will write a list of allowed \
gene IDs, allowed transcript IDs, and allowed gene and transcript IDs for mitochondrion.", required=True)
    parser.add_argument("--mito", help="Name of mitochondrial sequence. Default = chrM", default="chrM")
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
    
    # List of allowable biotypes from 10X Genomics:
    # https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps)

    allow_biotype = set(['protein_coding', 'lncRNA', 'IG_C_gene', 'IG_D_gene', \
        'IG_J_gene', 'IG_LV_gene', 'IG_V_gene', 'IG_V_pseudogene', 'IG_J_pseudogene', \
        'IG_C_pseudogene', 'TR_C_gene', 'TR_D_gene', 'TR_J_gene', 'TR_V_gene', \
        'TR_V_pseudogene', 'TR_J_pseudogene'])

    # List of non-allowable tags:
    
    disallow_tag = set(['readthrough_transcript', 'readthrough_gene', 'PAR', \
        'fragmented_locus', 'low_sequence_quality'])

    f_gz = False
    f = None
    if options.gtf[-3:] == '.gz':
        f_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f_gz = False
        f = open(options.gtf, 'r')
    
    out_gene = open('{}.genes'.format(options.output_prefix), 'w')
    out_tx = open('{}.tx'.format(options.output_prefix), 'w')
    out_gene_mito = open('{}.mito'.format(options.output_prefix), 'w')
    out_tx2gene = open('{}.tx2gene'.format(options.output_prefix), 'w')
    out_mito_gtf = open('{}.mito.gtf'.format(options.output_prefix), 'w')

    for line in f:
        if f_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != '#':
            dat = line.split('\t')
            
            if dat[2] == 'gene' or dat[2] == 'transcript':
                has_allowed_type = False
                has_disallowed_tag = False
                tags = {}
                for elt in dat[8].strip().rstrip(';').split(';'):
                    elt = elt.strip()
                    k, v = elt.split(' ')
                    v = v.strip('"')
                    tags[k] = v
                    if dat[2] == 'gene':
                        if k == 'gene_type' and v in allow_biotype:
                            has_allowed_type = True
                        elif k == 'tag' and v in disallow_tag:
                            has_disallowed_tag = True
                    elif dat[2] == 'transcript':
                        if k == 'transcript_type' and v in allow_biotype:
                            has_allowed_type = True
                        elif k == 'tag' and v in disallow_tag:
                            has_disallowed_tag = True
                
                if dat[2] == 'transcript':
                    print("{}\t{}".format(tags['transcript_id'].split('.')[0], \
                        tags['gene_id'].split('.')[0]), file=out_tx2gene)

                if has_allowed_type and not has_disallowed_tag:
                    if 'gene_id' in tags:
                        ensg = tags['gene_id'].split('.')[0]
                        name = tags['gene_name']
                        if name == ensg:
                            if ensg in ens2name:
                                name = ens2name[ensg]
                        elif ensg in ens2name and ens2name[ensg] != name:
                            name = ens2name[ensg]
                        
                        hgnc = ""
                        if 'hgnc_id' in tags:
                            hgnc = tags['hgnc_id']
                            if hgnc in hgnc2name and hgnc2name[hgnc] != name:
                                name = hgnc2name[hgnc]
                        
                        if dat[2] == 'gene':
                            if dat[0] == options.mito:
                                print("{}\t{}\t{}".format(ensg, name, hgnc), file=out_gene_mito)
                                # Print the actual annotation data to a file that can be lifted over
                                # separately (i.e. using liftOff)
                                tags['gene_name'] = name
                                dat[8] = join_tags(tags)
                                print("\t".join(dat), file=out_mito_gtf)
                            else:
                                print("{}\t{}\t{}".format(ensg, name, hgnc), file=out_gene)
                        else:
                            if 'transcript_id' in tags:
                                enst = tags['transcript_id'].split('.')[0]
                                if dat[0] != options.mito:
                                    print(enst, file=out_tx)
    f.close()
    out_gene.close()
    out_tx.close()
    out_gene_mito.close()
    out_tx2gene.close()
    out_mito_gtf.close()


if __name__ == '__main__':
    sys.exit(main(sys.argv))

