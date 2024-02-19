#! /usr/bin/env python3
import sys
import os
import argparse
import subprocess
"""
Program to run all of these other programs in a convenient way.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    
    # get_allowlist.py
    parser.add_argument("--gencode", help="Current Gencode annotation for human. \
This was tested on version 45. Can be gzipped.", required=True)
    parser.add_argument("--hgnc_ens", "-he", help="File with 3 (tab separated) columns: \
HGNC ID, HGNC approved symbol, Ensembl ID (or blank if missing). Required.", required=True)
    parser.add_argument("--output_prefix", "-o", help="Prefix for output files. The revised annotation \
will be [output_prefix].gtf.", required=True)
    parser.add_argument("--mito_name_hg38", help="Name of mitochondrial sequence in Gencode. Default = chrM", \
        default="chrM")
    
    # cleanup_cat_annotation.py
    parser.add_argument("--cat_gtf", help="CAT annotation to filter (GTF format). Can be gzipped.", required=True)
    parser.add_argument("--excl_contigs", "-ec", help="An optional list of contigs (target genome) to remove \
(should be a file)", required=False)
    parser.add_argument("--rename_contigs", "-rc", help="A file mapping old -> new contig (for target genome) \
names (tab separated, optional)", required=False)
    parser.add_argument("--drop_denovo", "-d", action="store_true", help="Set this option \
to remove all de novo predicted genes (not based on homology)")
    
    # Arguments to liftOff mitochondrial annotation
    mito_opts = parser.add_argument_group()
    mito_opts.add_argument("--cat_fasta", help="(bgzipped or not) FASTA file for species with CAT annotation (required)", required=True)
    mito_opts.add_argument("--hg38_fasta", help="(bgzipped or not) FASTA file for hg38 (required)", required=True)
    mito_opts.add_argument("--cat_mito", help="Name of mitochondrial sequence in CAT-annotated genome (required)", required=True)
    mito_opts.add_argument("--liftoff_path", help="Path to liftOff program (OPTIONAL; default = in path)", required=False)
    mito_opts.add_argument("--samtools_path", help="Path to samtools program (OPTIONAL; default = in path)", required=False)
    mito_opts.add_argument("--minimap2_path", help="Path to minimap2 program (OPTIONAL; default = in path)", required=False)

    # Arguments for rescuing genes from an old Ensembl annotation (optional)
    rescue_opts = parser.add_argument_group()
    rescue_opts.add_argument("--lift_old_annotation", '-l', help="If genes are eliminated from the CAT annotation \
because all their transcripts have been eliminated, but the genes themselves pass filter, it may be possible to \
rescue annotations for those genes by lifting over an Ensembl annotation from a previous assembly of that species \
to the current assembly. This requires an Ensembl annotation for the species, a UCSC chain file mapping between the \
old and new version of the assembly, and the programs gtfToGenePred, liftOver, and genePredToGtf.", required=False, \
        action="store_true")
    rescue_opts.add_argument("--species_id", "-s", help="Species identifier for CAT-annotated genome in Ensembl \
databases. This shoud be lowercase first letter of species plus full lowercase genus, with no spaces. For example, \
chimpanzee is ptroglodytes, and rhesus macaque is mmulatta. (required if -l set)", required=False)
    rescue_opts.add_argument("--old_ens_gtf", help="Ensembl annotation on an old assembly of the same species for which \
the CAT annotation exists (required if -l set). Can be gzipped, and should be the \"chr\" annotation (i.e. \
the one that does not include alt haplotypes and unplaced scaffolds).", required=False)
    rescue_opts.add_argument("--chain", help="UCSC chain file mapping from old assembly of this species to the current \
one (required if -l set)", required=False)
    rescue_opts.add_argument("--gtfToGenePred_path", help="Path to gtfToGenePred (OPTIONAL; default = in path)", required=False)
    rescue_opts.add_argument("--liftOver_path", help="Path to liftOver (OPTIONAL; default = in path)", required=False)
    rescue_opts.add_argument("--genePredToGtf_path", help="Path to genePredToGtf (OPTIONAL; default = in path)", required=False)
    
    parsed = parser.parse_args()
    if parsed.lift_old_annotation:
        if parsed.species_id is None:
            print("ERROR: --species_id is required if lifting an old annotation", file=sys.stderr)
            exit(1)
        if parsed.old_ens_gtf is None:
            print("ERROR: --old_ens_gtf is required if lifting an old annotation", file=sys.stderr)
            exit(1)
        if parsed.chain is None:
            print("ERROR: --chain is required if lifting an old annotation", file=sys.stderr)
            exit(1)
    return parser.parse_args()

def main(args):
    options = parse_args()
    curdir = os.path.dirname(os.path.realpath(__file__))
    
    # Get list of allowed genes
    cmd1 = ['{}/get_allowlist.py'.format(curdir), '--gtf', options.gencode, '-he', options.hgnc_ens, \
        '-o', options.output_prefix, '--mito', options.mito_name_hg38]
    print(" ".join(cmd1), file=sys.stderr)
    subprocess.call(cmd1)
    
    # Filter the CAT annotation and remove mitochondrial genes
    cmd2 = ['{}/cleanup_cat_annotation.py'.format(curdir), '--gtf', options.cat_gtf, '-a', options.output_prefix, \
        '-o', '{}'.format(options.output_prefix)]
    if options.excl_contigs is not None:
        cmd2.append('--excl_contigs')
        cmd2.append(options.excl_contigs)
    if options.rename_contigs is not None:
        cmd2.append('--rename_contigs')
        cmd2.append(options.rename_contigs)
    if options.drop_denovo:
        cmd2.append("--drop_denovo")
    print(" ".join(cmd2), file=sys.stderr)
    subprocess.call(cmd2)

    # Use liftoff to create new mitochondrial annotation
    cmd3 = ['{}/mito_liftoff.py'.format(curdir), '--cat_fasta', options.cat_fasta, '--hg38_fasta', options.hg38_fasta, \
        '--cat_mito', options.cat_mito, '--allowlist_base', options.output_prefix, '-o', options.output_prefix]
    if options.liftoff_path is not None:
        cmd3.append('--liftoff_path')
        cmd3.append(options.liftoff_path)
    if options.samtools_path is not None:
        cmd3.append('--samtools_path')
        cmd3.append(options.samtools_path)
    if options.minimap2_path is not None:
        cmd3.append("--minimap2_path")
        cmd3.append(options.minimap2_path)
    print(" ".join(cmd3), file=sys.stderr)
    subprocess.call(cmd3)
    
    # Clean up
    os.unlink('{}.mito'.format(options.output_prefix))
    os.unlink('{}.mito.gtf'.format(options.output_prefix))

    if options.species_id is not None:
        # Attempt to rescue by lifting from an old Ensembl annotation
        cmd4 = ['{}/get_homologous_genelist.R'.format(curdir), '{}.rm.genes.notx'.format(options.output_prefix), \
            options.species_id]
        f = open('{}.homology'.format(options.output_prefix), 'w')
        print(" ".join(cmd4), file=sys.stderr)
        p = subprocess.Popen(cmd4, stdout=f)
        out, err = p.communicate()
        f.close()
        cmd5 = ['{}/lift_old_ensembl_annotation.py'.format(curdir), '-G', '{}.homology'.format(options.output_prefix), \
            '-g', options.old_ens_gtf, '-c', options.chain, '--hgnc_ens', options.hgnc_ens, '-o', options.output_prefix]
        if options.gtfToGenePred_path is not None:
            cmd5.append("--gtfToGenePred_path")
            cmd5.append(options.gtfToGenePred_path)
        if options.liftOver_path is not None:
            cmd5.append("--liftOver_path")
            cmd5.append(options.liftOver_path)
        if options.genePredToGtf_path is not None:
            cmd5.append("--genePredToGtf_path")
            cmd5.append(options.genePredToGtf_path)

        print(" ".join(cmd5), file=sys.stderr)
        subprocess.call(cmd5)
        
        os.unlink('{}.homology'.format(options.output_prefix))

    # Clean up
    if os.path.isfile('{}.tx'.format(options.output_prefix)):
        os.unlink('{}.tx'.format(options.output_prefix))
    if os.path.isfile('{}.tx2gene'.format(options.output_prefix)):
        os.unlink('{}.tx2gene'.format(options.output_prefix))
    if os.path.isfile('{}.genes'.format(options.output_prefix)):
        os.unlink('{}.genes'.format(options.output_prefix))

    # Now put everything together into one annotation.
    


if __name__ == '__main__':
    sys.exit(main(sys.argv))
