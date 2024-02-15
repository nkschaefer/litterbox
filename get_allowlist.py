#! /usr/bin/env python3
import sys
import os
import argparse
import gzip
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
    parser.add_argument("--mito", "-m", help="Set to print a list of mitochondrial gene \
IDs. Default: print an allow list of other (non-mitochondrial) genes and their preferred \
names.", action="store_true")
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
    
    disallow_tag = set(['readthrough_transcript', 'PAR'])

    mito = set(['chrM', 'chrMT', 'MT', 'M'])

    f_gz = False
    f = None
    if options.gtf[-3:] == '.gz':
        f_gz = True
        f = gzip.open(options.gtf, 'r')
    else:
        f_gz = False
        f = open(options.gtf, 'r')

    for line in f:
        if f_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != '#':
            dat = line.split('\t')
            
            if (options.mito and dat[0] in mito) or \
                (not options.mito and dat[0] not in mito):
                if dat[2] == 'gene':
                    has_allowed_type = False
                    has_disallowed_tag = False
                    tags = {}
                    for elt in dat[8].strip().rstrip(';').split(';'):
                        elt = elt.strip()
                        k, v = elt.split(' ')
                        v = v.strip('"')
                        tags[k] = v
                        if k == 'gene_type' and v in allow_biotype:
                            has_allowed_type = True
                        elif k == 'tag' and v in disallow_tag:
                            has_disallowed_tag = True

                    if has_allowed_type and not has_disallowed_tag:
                        if 'gene_id' in tags:
                            ensg = tags['gene_id'].split('.')[0]
                            name = tags['gene_name']
                            if name == ensg:
                                if ensg in ens2name:
                                    name = ens2name[ensg]
                            elif ensg in ens2name and ens2name[ensg] != name:
                                name = ens2name[ensg]
                            
                            havana = ""
                            if 'havana_gene' in tags:
                                havana = tags['havana_gene'].split('.')[0]
                            hgnc = ""
                            if 'hgnc_id' in tags:
                                hgnc = tags['hgnc_id']
                                if hgnc in hgnc2name and hgnc2name[hgnc] != name:
                                    name = hgnc2name[hgnc]

                            print("{}\t{}\t{}\t{}".format(ensg, name, hgnc, havana))
    
    f.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))

