#! /usr/bin/env python3
import gffutils
import sys
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Convert GFF3 to GTF, preserving order so cellranger-arc mkref doesn't complain.")
    parser.add_argument("--gff3", "-g", help="Input GFF3 file", required=True)
    parser.add_argument("--output_prefix", "-o", help="Prefix for output files", required=True)
    return parser.parse_args()


gtf_dialect = {'field separator': '; ', \
    'fmt': 'gtf', \
    'keyval separator': ' ', \
    'leading semicolon': False, \
    'multival separator': ',', \
    'quoted GFF2 values': True, \
    'repeated keys': False, \
    'trailing semicolon': True}


def recurse_children(db, child):
    
    if child.featuretype == "intron":
        return

    #if "Parent" in child.attributes:
    #    del child.attributes["Parent"]
    """
    if child.featuretype != 'gene':
        print(child.featuretype)
        print(child.attributes)
        exit()
    """
    child.dialect = gtf_dialect
    
    # This stuff is needed by cellranger-atac mkref
    if 'gene_name' not in child.attributes:
        if 'Name' in child.attributes:
            child.attributes['gene_name'] = child.attributes['Name']
    if 'gene_type' not in child.attributes:
        if 'gene_biotype' in child.attributes:
            child.attributes['gene_type'] = child.attributes['gene_biotype']
    if 'gene_biotype' in child.attributes:
        del child.attributes['gene_biotype']

    print(child)
    recurse_children.printed.add(child.attributes['ID'][0])

    for child2 in db.children(child, level=1):
        if not child2.attributes["ID"][0] in recurse_children.printed:
            recurse_children(db, child2)

recurse_children.printed = set([])

def main(args):
    options = parse_args()
    fnbase = options.gff3.split('/')[-1]
    fnbase = fnbase.split('.')[0]
    db = None
    if os.path.isfile('{}.db'.format(options.output_prefix)):
        db = gffutils.interface.FeatureDB('{}.db'.format(options.output_prefix))
    else:
        db = gffutils.create_db(options.gff3, dbfn="{}.db".format(options.output_prefix), force=True, \
            keep_order=False, sort_attribute_values=False, merge_strategy="create_unique")
    for gene in db.features_of_type('gene'):
        recurse_children(db, gene)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
