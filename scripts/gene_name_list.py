#! /usr/bin/env python3
import sys
import os
from cleanup_cat_annotation import get_tags
"""
Get a list of all gene IDs and names, tab separated, from a GTF.
Can be used to check if any names map to > 1 gene.
"""
def main(args):
    for line in sys.stdin:
        line = line.rstrip()
        dat = line.split('\t')
        if dat[2] == 'gene':
            tags = get_tags(dat)
            if 'gene_name' in tags:
                name = tags['gene_name']
                if 'gene_id' in tags:
                    gid = tags['gene_id']
                    print("{}\t{}".format(name, gid))

if __name__ == '__main__':
    sys.exit(main(sys.argv))
