#! /usr/bin/env python3
import os
import sys
import argparse
import subprocess
"""
Given mitochondrial annotations (from Gencode), and mitochondrial genome sequences
from hg38 (Gencode target) and the species of interest, uses liftOff to map the
mitochondrial annotations onto the genome of the other species.
"""
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cat_fasta", 
        help="(bgzipped or not) FASTA file for species with CAT annotation", 
        required=True
    )
    parser.add_argument(
        "--hg38_fasta",
        help="(bgzipped or not) FASTA file for hg38",
        required=True
    )
    parser.add_argument(
        "--liftoff_path", 
        help="Path to liftOff program (OPTIONAL; default = in path)",
        required=False
    )
    parser.add_argument(
        "--samtools_path",
        help="Path to samtools program (OPTIONAL; default = in path)",
        required=False
    )
    parser.add_argument(
        "--minimap2_path",
        help="Path to minimap2 program (OPTIONAL; default = in path)",
        required=False
    )
    parser.add_argument(
        "--allowlist_base", 
        '-a', 
        help="Output prefix used for get_allowlist.py", 
        required=True
    )
    parser.add_argument(
        "--hg38_mito", 
        help="Name of mitochondrion in hg38", 
        required=False, 
        default="chrM"
    )
    parser.add_argument(
        "--cat_mito", 
        help="Name of mitochondrion in CAT", 
        required=False, 
        default=None
    )
    parser.add_argument(
        "--output_prefix", 
        "-o", 
        help="Prefix for output files to create", 
        required=True
    )

    return parser.parse_args()

def main(args):
    options = parse_args()
    
    samtools = 'samtools'
    if options.samtools_path is not None:
        if options.samtools_path[-1] == '/':
            samtools = options.samtools_path + 'samtools'
        else:
            samtools = options.samtools_path + '/samtools'

    liftoff = 'liftoff'
    if options.liftoff_path is not None:
        if options.liftoff_path[-1] == '/':
            liftoff = options.liftoff_path + 'liftoff'
        else:
            liftoff = options.liftoff_path + '/liftoff'

    mm2 = None
    if options.minimap2_path is not None:
        if options.minimap2_path[-1] == '/':
            mm2 = options.minimap2_path = 'minimap2'
        else:
            mm2 = options.minimap2_path + "/minimap2"

    # Make sure both genomes are samtools faidxed
    if not os.path.isfile('{}.fai'.format(options.hg38_fasta)):
        print("Indexing {}...".format(options.hg38_fasta), file=sys.stderr)
        subprocess.call([samtools, 'faidx', options.hg38_fasta])
    if not os.path.isfile('{}.fai'.format(options.cat_fasta)):
        print("Indexing {}...".format(options.cat_fasta), file=sys.stderr)
        subprocess.call([samtools, 'faidx', options.cat_fasta])

    # Make sure we know the name of the mitochondrion
    cat_mito = None
    if options.cat_mito is None:
        f = open('{}.fai'.format(options.cat_fasta), 'r')
        seqs = set([])
        for line in f:
            line = line.rstrip()
            dat = line.split('\t')
            seqs.add(dat[0])
        f.close()
        candidates = ['chrM', 'M', 'MT', 'chrMT']
        for cand in candidates:
            if cand in seqs:
                print("Using {} as mitochondrial seq".format(cand), \
                    file=sys.stderr)
                cat_mito = cand
                break
        if cat_mito is None:
            print("ERROR: could not guess mitochondrial sequence name in {}".\
                format(options.cat_fasta), file=sys.stderr)
            print("Please provide manually.", file=sys.stderr)
            exit(1)
    else:
        cat_mito = options.cat_mito

    # Get mitochondrial genomes
    f = open('{}.mito.gencode.fa'.format(options.output_prefix), 'w')
    p = subprocess.Popen([
        samtools, 
        'faidx', 
        options.hg38_fasta, 
        options.hg38_mito
    ], stdout=f)
    out, err = p.communicate()
    f.close()
    
    f = open('{}.mito.cat.fa'.format(options.output_prefix), 'w')
    p = subprocess.Popen([
        samtools, 
        'faidx', 
        options.cat_fasta, 
        cat_mito], 
    stdout=f)
    out, err = p.communicate()
    f.close()

    # Run liftOff
    cmd = [
        liftoff, 
        '-g', '{}.mito.gtf'.format(options.allowlist_base), 
        '-o', '{}.mito.liftoff.gtf'.format(options.output_prefix),
        '-u', '{}.mito.unmapped'.format(options.output_prefix),
        '-flank', '10'
    ]
    
    if mm2 is not None:
        cmd.append('-m')
        cmd.append(mm2)
    
    cmd.append('{}.mito.cat.fa'.format(options.output_prefix))
    cmd.append('{}.mito.gencode.fa'.format(options.output_prefix))
    print(" ".join(cmd), file=sys.stderr)
    subprocess.call(cmd)
    
if __name__ == '__main__':
    sys.exit(main(sys.argv))
