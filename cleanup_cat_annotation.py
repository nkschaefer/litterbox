#! /usr/bin/env python3
import sys
import os
from collections import Counter
import argparse
import gzip
"""
After converting from GFF3 to GTF, filter the GTF by removing any genes 
based on human annotations that would be removed according to 10X Genomics's
guidelines. Also remove any de novo gene predictions that are not likely to 
be accurate/useful (poor_alignment and possible_paralog with no RNA support).

Also removes mitochondrial annotations (unreliable) and all mitochondrial genes
annotated to the be on the autosomes (unreliable or NUMTs, which will soak up
true mitochondrial reads as well).

Optionally can exclude specific contigs (i.e. random/unplaced scaffolds) and/or
rename contigs.

After this, still need to add in mitochondrial annotations.
"""

def parse_args():
    parser = argparse.ArgumentParser(description="__doc__")
    parser.add_argument("--gtf", "-g", help="GTF to filter", required=True)
    parser.add_argument("--allowlist", "-a", help="Filtered list of allowed ENSG and name, \
tab separated, plus optional additional columns at the end (ignored)", required=True)
    parser.add_argument("--mito_genes", "-m", help="A list of mitochondrial gene IDs to \
remove (ENSG IDs).", required=True)
    parser.add_argument("--excl_contigs", "-ec", help="An optional list of contigs to remove \
(should be a file)", required=False)
    parser.add_argument("--rename_contigs", "-rc", help="A file mapping old -> new contig \
names (tab separated, optional)", required=False)
    return parser.parse_args()

def get_tags(dat):
    tags = {}
    for elt in dat[8].strip().rstrip(';').split(';'):
        elt = elt.strip()
        k, v = elt.split(' ')
        v = v.strip('"')
        tags[k] = v
    return tags

def join_tags(tags):
    tagprint = []
    for k in tags:
        tagprint.append('{} "{}"'.format(k, tags[k]))
    return "; ".join(tagprint) + ";"

def main(args):
    
    options = parse_args()

    g2n = {}
    f = open(options.allowlist, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        ensg = dat[0]
        name = dat[1]
        g2n[ensg] = name
    f.close()

    mito = set([])
    f = open(options.mito_genes, 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        mito.add(dat[0])
    f.close()

    contigs_excl = set([])
    if options.excl_contigs is not None:
        f = open(options.excl_contigs, 'r')
        for line in f:
            line = line.rstrip()
            contigs_excl.add(line)
        f.close()
    
    contigs_rename = {}
    if options.rename_contigs is not None:
        f = open(options.rename_contigs, 'r')
        for line in f:
            line = line.rstrip()
            if len(line) > 0:
                dat = line.split('\t')
                contigs_rename[dat[0]] = dat[1]
        f.close()

    # Also want to exclude mitochondrial sequences, since we have to re-do them
    contigs_excl.add('chrM')
    contigs_excl.add('MT')
    contigs_excl.add('M')
    contigs_excl.add('chrMT')

    gene_txcounts = Counter()
    tx_rm = set([])

    tx_disallowed = set(['possible_paralog', 'poor_alignment'])

    gene_cand_name = {}
    
    # Make two passes
    f = None
    f_gz = False
    if options.gtf[-3:] == '.gz':
        f = gzip.open(options.gtf, 'r')
        f_gz = True
    else:
        f = open(options.gtf, 'r')
        f_gz = False

    for line in f:
        if f_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            
            if dat[0] in contigs_excl:
                continue
            tags = get_tags(dat)            
            if dat[2] == 'gene':
                if 'source_gene' in tags and tags['source_gene'] != 'None':
                    # Look up ENS ID to see if it passes
                    ensg = tags['source_gene']
                    ensg = ensg.split('.')[0]
                    
                    if ensg not in mito and ensg in g2n:
                        gene_cand_name[tags['gene_id']] = g2n[ensg]
                else:
                    # Store name for gene and wait
                    gene_cand_name[tags['gene_id']] = tags['gene_name']
            
            elif dat[2] == 'transcript':
                if 'source_gene' not in tags or tags['source_gene'] == 'None':
                    # This is predicted transcript not based on homology
                    txrm = False
                    if 'transcript_class' not in tags:
                        txrm = True
                    elif tags['transcript_class'] == 'poor_alignment':
                        txrm = True
                    elif tags['transcript_class'] == 'possible_paralog':
                        # Check the exon RNA support field.
                        if 'exon_rna_support' not in tags:
                            txrm = True
                        else:
                            ref_support = False
                            if 'reference_support' in tags and tags['reference_support'] \
                                == 'True':
                                ref_support = True
                            rna_support = False
                            if 'rna_support' in tags and tags['rna_support'] == 'True':
                                rna_support = True
                            pacbio_support = False
                            if 'pacbio_isoform_supported' in tags and \
                                tags['pacbio_isoform_supported'] == 'True':
                                pacbio_support = True

                            if not rna_support and not pacbio_support:
                                txrm = True
                    if txrm:
                        tx_rm.add(tags['transcript_id'])
                    else:
                        gid = tags['Parent']
                        gene_txcounts[gid] += 1
                else:
                    # Gene will be removed or not based on Ensembl annotations
                    gid = tags['Parent']
                    if 'source_gene' not in tags or tags['source_gene'].split('.')[0] not in mito: 
                        gene_txcounts[gid] += 1
    f.close()

    # Now, do a second pass and only keep genes with at least one valid transcript.
    # Also only keep CDSs and exons whose parent is a valid transcript.
    if f_gz:
        f = gzip.open(options.gtf, 'r')
    else:
        f = open(options.gtf, 'r')
    for line in f:
        if f_gz:
            line = line.decode().rstrip()
        else:
            line = line.rstrip()
        dat = line.split('\t')

        if dat[0] in contigs_excl:
            continue

        tags = get_tags(dat)
        
        printLine = False
        if dat[2] == 'gene':
            gid = tags['gene_id']
            if gid in gene_cand_name and gene_txcounts[gid] > 0:
                # Pass
                tags['gene_name'] = gene_cand_name[gid]
                printLine = True
        elif dat[2] == 'transcript':
            # Gene must pass AND transcript must pass
            gid = tags['gene_id']
            if gid in gene_cand_name and gene_txcounts[gid] > 0:
                if tags['transcript_id'] not in tx_rm:
                    printLine = True
        else:
            # Everything else is a child of transcript.
            gid = tags['gene_id']
            tid = tags['transcript_id']
            if gid in gene_cand_name and gene_txcounts[gid] > 0 and tid not in tx_rm:
                printLine = True
        if printLine:
            if dat[0] in contigs_rename:
                dat[0] = contigs_rename[dat[0]]
            # Stick in gene name to use
            tags['gene_name'] = gene_cand_name[tags['gene_id']]
            tags['Name'] = gene_cand_name[tags['gene_id']]
            dat[8] = join_tags(tags)
            print("\t".join(dat))


if __name__ == '__main__':
    sys.exit(main(sys.argv))
