#! /usr/bin/env python3
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", "-g", help="GTF from liftoff to filter", required=True)
    parser.add_argument("--chrom", "-c", help="True chromosome for alignment", required=True)
    parser.add_argument("--start_pos", "-p", help="Start position of genomic segment aligned to", type=int,
            required=True)
    parser.add_argument("--out_gtf", "-o", help="Name of output GTF to create", required=True)
    parser.add_argument("--out_genes", "-G", help="Name of list of rescued genes to create", required=True)
    parser.add_argument("--gene_id", '-i', help="Gene ID field in GTF", default="gene_id", required=False)
    parser.add_argument("--gene_name", '-n', help="Gene name field in GTF", default="gene_name", required=False)
    return parser.parse_args()

def main(args):
    options = parse_args()
    f_genes = open(options.out_genes, 'w')
    f_gtf = open(options.out_gtf, 'w')
    f = open(options.gtf, 'r')
    genes_rescue = set([])
    for line in f:
        line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            dat[0] = options.chrom
            dat[3] = str(int(dat[3]) + options.start_pos)
            dat[4] = str(int(dat[4]) + options.start_pos)
            gid = None
            gname = None
            for elt in dat[8].split('; '):
                elt = elt.strip(';')
                if len(elt) > 0:
                    eltsplit = elt.split(' ')
                    if len(eltsplit) == 2:
                        k, v = eltsplit
                        v = v.strip('"')
                        if k == options.gene_id:
                            gid = v
                        elif k == options.gene_name:
                            gname = v
            if gname is not None:
                genes_rescue.add(gname)
            print("\t".join(dat), file=f_gtf)
    f.close()
    f_gtf.close()
    for g in sorted(list(genes_rescue)):
        print(g, file=f_genes)
    f_genes.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
