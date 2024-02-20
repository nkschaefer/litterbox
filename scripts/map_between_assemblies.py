#! /usr/bin/env python3
import sys
import os
"""
The HAL annotations sometimes have weird, different names
for chromosomes/scaffolds/contigs than the official assemblies.

This script takes two fasta index files, one for each version
of the same assembly, and matches up sequences by length, then
outputs a map of sequence in fai 1 to sequence in fai 2.
"""

def main(args):
    if len(args) < 3:
        print("USAGE: map_between_assemblies.py file1.fai file2.fai", file=sys.stderr)
        exit(1)
    len2seq1 = {}
    f = open(args[1], 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        len2seq1[int(dat[1])] = dat[0]
    f.close()

    f = open(args[2], 'r')
    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        length = int(dat[1])
        if length in len2seq1:
            print("{}\t{}".format(len2seq1[length], dat[0]))
        else:
            print("WARNING: {} {} does not match {}".format(dat[0], length, args[1]),
                file=sys.stderr)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
