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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gtf", 
        "-g", 
        help="GTF to filter", 
        required=True
    )
    parser.add_argument(
        "--allowlist_base", 
        "-a", 
        help="Output prefix from get_allowlist.py run",
        required=True
    )
    parser.add_argument(
        "--excl_contigs",
        "-ec", 
        help="An optional list of contigs to remove (should be a file)", 
        required=False
    )
    parser.add_argument(
        "--rename_contigs", 
        "-rc", 
        help="A file mapping old -> new contig names (tab separated, optional)", 
        required=False
    )
    parser.add_argument(
        "--drop_denovo", 
        "-d", 
        action="store_true", 
        help="Set this option to remove all de novo predicted genes (not based \
on homology)"
    )
    parser.add_argument(
        "--output_prefix", 
        "-o", 
        required=True, 
        help="Prefix for output files"
    )
    return parser.parse_args()

def get_tags(dat):
    tags = {}
    try:
        for elt in dat[8].strip().rstrip(';').split(';'):
            elt = elt.strip()
            k, v = elt.split(' ')
            v = v.strip('"')
            tags[k] = v
    except e:
        print("ERROR: input likely not GTF format.", file=sys.stderr)
        exit(1)
    return tags

def join_tags(tags):
    tagprint = []
    for k in tags:
        tagprint.append('{} "{}"'.format(k, tags[k]))
    return "; ".join(tagprint) + ";"

def main(args):
    
    options = parse_args()

    g2n = {}
    f = open('{}.genes'.format(options.allowlist_base), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        ensg = dat[0]
        name = dat[1]
        g2n[ensg] = name
    f.close()

    tx_pass = set([])
    f = open('{}.tx'.format(options.allowlist_base), 'r')
    for line in f:
        line = line.rstrip()
        tx_pass.add(line)
    f.close()

    mito = set([])
    f = open('{}.mito'.format(options.allowlist_base), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        mito.add(dat[0])
    f.close()
    
    tx2gene = {}
    f = open('{}.tx2gene'.format(options.allowlist_base), 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        tx2gene[dat[0]] = dat[1]
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

    gene_txcounts = Counter()
    tx_rm = set([])

    gene_cand_name = {}
    
    gene2source = {}

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
                    gene2source[tags['gene_id']] = ensg
                    
                    if ensg not in mito and ensg in g2n:
                        gene_cand_name[tags['gene_id']] = g2n[ensg]
                
                else:
                    # Store name for gene and wait
                    if 'gene_name' in tags:
                        gene_cand_name[tags['gene_id']] = tags['gene_name']
                    elif 'source_gene_common_name' in tags:
                        gene_cand_name[tags['gene_id']] = \
                            tags['source_gene_common_name'].split('.')[0]
                    else:
                        gene_cand_name[tags['gene_id']] = tags['gene_id']

            elif dat[2] == 'transcript':
                txrm = False
                if 'source_transcript' in tags:
                    enst = tags['source_transcript'].split('.')[0]
                    if enst not in tx_pass:
                        txrm = True
                elif 'source_gene' not in tags or tags['source_gene'] == 'None':
                    # This is predicted transcript not based on homology
                    if 'transcript_class' not in tags:
                        txrm = True
                    elif tags['transcript_class'] == 'poor_alignment':
                        txrm = True
                    elif options.drop_denovo and tags['transcript_class'] in \
                        ['possible_paralog', 'putative_novel_isoform', 'putative_novel']:
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
                            
                            if ref_support and (rna_support or pacbio_support):
                                # Keep it
                                pass
                            else:
                                txrm = True
                if txrm:
                    tx_rm.add(tags['transcript_id'])
                else:
                    gid = tags['gene_id']
                    gene_txcounts[gid] += 1
    f.close()

    out_gtf = open("{}.main.gtf".format(options.output_prefix), 'w')
    out_gene_drop = open("{}.rm.genes".format(options.output_prefix), 'w')
    out_gene_drop_tx = open("{}.rm.genes.notx".format(options.output_prefix), 'w')
    out_tx_drop = open("{}.rm.tx".format(options.output_prefix), 'w')

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

        if tags['gene_id'] in gene2source and gene2source[tags['gene_id']] != 'nan':
            tags['source_gene'] = gene2source[tags['gene_id']]
        
        if dat[2] == 'gene':
            gid = tags['gene_id']
            rm_because_tx = False
            if gid in gene_cand_name and gene_txcounts[gid] > 0:
                # Pass
                tags['gene_name'] = gene_cand_name[gid]
                printLine = True
            else:
                if gid in gene_cand_name:
                    rm_because_tx = True
                name = gid
                if gid in gene_cand_name:
                    name = gene_cand_name[gid]
                elif gid in gene2source:
                    eg = gene2source[gid]
                    if eg in g2n:
                        name = g2n[eg]
                source = ""
                if 'source_gene' in tags and tags['source_gene'] != 'nan':
                    source = tags['source_gene']
                print("{}\t{}\t{}".format(gid, name, source), file=out_gene_drop)
                if rm_because_tx:
                    print("{}\t{}\t{}".format(gid, name, source), file=out_gene_drop_tx)

        elif dat[2] == 'transcript':
            # Gene must pass AND transcript must pass
            gid = tags['gene_id']
            if gid in gene_cand_name and gene_txcounts[gid] > 0:
                if tags['transcript_id'] not in tx_rm:
                    printLine = True
                else:
                    print(tags['transcript_id'], file=out_tx_drop)
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
            dat[8] = join_tags(tags)
            print("\t".join(dat), file=out_gtf)

    out_gtf.close()
    out_gene_drop.close()
    out_gene_drop_tx.close()
    out_tx_drop.close()



if __name__ == '__main__':
    sys.exit(main(sys.argv))
