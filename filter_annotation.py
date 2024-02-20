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
    parser.add_argument(
        "--gencode",
        help="Current Gencode annotation for human. This was tested on version \
45. Can be gzipped.", 
        required=True
    )
    parser.add_argument(
        "--hgnc_ens",
        "-he", 
        help="File with 3 (tab separated) columns: HGNC ID, HGNC approved symbol, \
Ensembl ID (or blank if missing). OPTIONAL; default = data/hgncid_name_ens.txt.", 
        required=False,
    )
    parser.add_argument(
        "--output_prefix", 
        "-o", 
        help="Prefix for output files. The revised annotation will be \
[output_prefix].gtf.", 
        required=True
    )
    parser.add_argument(
        "--mito_name_hg38", 
        help="Name of mitochondrial sequence in Gencode. Default = chrM",
        default="chrM"
    )
    
    # filter_human_annotation.py
    parser.add_argument(
        "--human_annotation",
        help="Optionally filter a human annotation (e.g. 10X) to only include genes \
and transcripts in the filtered list, and to use the HGNC-preferred names.",
        required=False
    )

    # cleanup_cat_annotation.py
    parser.add_argument(
        "--cat_gff3", 
        help="CAT annotation to filter (GFF3 format). Can be gzipped.", 
        required=True
    )
    parser.add_argument(
        "--excl_contigs", 
        "-ec", 
        help="An optional list of contigs (target genome) to remove (should be a file)", 
        required=False
    )
    parser.add_argument(
        "--rename_contigs", 
        "-rc", 
        help="A file mapping old -> new contig (for target genome) names (tab \
separated, optional)", 
        required=False
    )
    parser.add_argument(
        "--drop_denovo", 
        "-d", 
        action="store_true", 
        help="Set this option to remove all de novo predicted genes (not \
based on homology)"
    )
    
    # Arguments to liftoff mitochondrial annotation
    mito_opts = parser.add_argument_group()
    mito_opts.add_argument(
        "--cat_fasta", 
        help="(bgzipped or not) FASTA file for species with CAT annotation (required)", 
        required=True
    )
    mito_opts.add_argument(
        "--hg38_fasta", 
        help="(bgzipped or not) FASTA file for hg38 (required)", 
        required=True
    )
    mito_opts.add_argument(
        "--cat_mito", 
        help="Name of mitochondrial sequence in CAT-annotated genome (required)", 
        required=True
    )
    mito_opts.add_argument(
        "--liftoff_path", 
        help="Path to liftOff program (OPTIONAL; default = in path)", 
        required=False
    )
    mito_opts.add_argument(
        "--samtools_path", 
        help="Path to samtools program (OPTIONAL; default = in path)", 
        required=False
    )
    mito_opts.add_argument(
        "--minimap2_path", 
        help="Path to minimap2 program (OPTIONAL; default = in path)", 
        required=False
    )

    # Arguments for rescuing genes from an old Ensembl annotation (optional)
    rescue_opts = parser.add_argument_group()
    rescue_opts.add_argument(
        "--scavenge_ensembl_annotation", 
        help="If genes are eliminated from the CAT annotation because all \
their transcripts have been eliminated, but the genes themselves pass filter,\
it may be possible to rescue annotations for those genes by pulling them from \
an Ensembl annotation on the same assembly, or lifting over an Ensembl \
annotation from a previous assembly of that species to the current assembly. \
This requires an Ensembl annotation for the species, and if the assemblies \
are different then it also requires  a UCSC chain file mapping between the \
old and new version of the assembly, and the programs gtfToGenePred, liftOver, \
and genePredToGtf.", 
        required=False, 
        action="store_true"
    )
    rescue_opts.add_argument(
        "--species_id", 
        "-s", 
        help="Species identifier for CAT-annotated genome in Ensembl databases. \
This shoud be lowercase first letter of species plus full lowercase genus, with \
no spaces. For example, chimpanzee is ptroglodytes, and rhesus macaque is mmulatta. \
(required if -l set)", 
        required=False
    )
    rescue_opts.add_argument(
        "--old_ens_gtf", 
        help="Ensembl annotation on an old assembly of the same species for which \
the CAT annotation exists (required if -l set). Can be gzipped, and should be the \
\"chr\" annotation (i.e. the one that excludes alt haplotypes and unplaced scaffolds).", 
        required=False
    )
    rescue_opts.add_argument(
        "--ensembl_fasta",
        help="FASTA file for Ensembl genome, so we can alter sequence names as \
necessary. Required if --scavenge_ensembl_annotation and the Ensembl assembly is \
the same version as the one used for CAT (i.e. you are not providing a --chain file)",
        required=False
    )
    rescue_opts.add_argument(
        "--ucsc_fasta",
        help="FASTA file for UCSC version of the genome the CAT annotation corresponds \
to. The HAL-aligned genomes the CAT annotation used often have renamed chromosomes, \
and this will ensure that sequence names come out matching the CAT annotation. Required \
only if the Ensembl annotation is on a different assembly and you are lifting over to \
the assembly CAT used (i.e. you are providing a --chain file)",
        required=False
    )
    rescue_opts.add_argument(
        "--chain", 
        help="UCSC chain file mapping from old assembly of this species to the current \
one (required if --scavenge_ensembl_annotation and Ensembl assembly is different version \
than the one used for CAT)",
        required=False
    )
    rescue_opts.add_argument(
        "--gtfToGenePred_path",
        help="Path to gtfToGenePred (OPTIONAL; default = in path)",
        required=False
    )
    rescue_opts.add_argument(
        "--liftOver_path", 
        help="Path to liftOver (OPTIONAL; default = in path)",
        required=False
    )
    rescue_opts.add_argument(
        "--genePredToGtf_path",
        help="Path to genePredToGtf (OPTIONAL; default = in path)",
        required=False
    )
    
    parsed = parser.parse_args()
    if parsed.scavenge_ensembl_annotation:
        if parsed.species_id is None:
            print(
                "ERROR: --species_id is required if lifting an old annotation", 
                file=sys.stderr
            )
            exit(1)
        if parsed.old_ens_gtf is None:
            print(
                "ERROR: --old_ens_gtf is required if lifting an old annotation",
                file=sys.stderr
            )
            exit(1)
        if parsed.chain is None:
            print(
                "WARNING: without chain file (--chain), the Ensembl annotation \
{} MUST be on the same genome assembly as the CAT annotation.".format(parsed.old_ens_gtf),
                file=sys.stderr
            )
        if parsed.chain is None and parsed.ensembl_fasta is None:
            print(
                "ERROR: you must provide an --ensembl_fasta genome if scavenging \
an annotation on the same assembly", file=sys.stderr
            )
        if parsed.chain is not None and parsed.ucsc_fasta is None:
            print(
                "ERROR: you must provide a --ucsc_fasta genome if scavenging \
an Ensembl annotation on a different assembly, then lifting over to a UCSC versioned \
assembly",
                file=sys.stderr
            )
    return parser.parse_args()

def main(args):
    options = parse_args()
    curdir = os.path.dirname(os.path.realpath(__file__))
    
    # If the user didn't specify a HGNC ID -> Approved symbol -> Ensembl ID
    # mapping file, use the one packaged here
    hgnc_ens = options.hgnc_ens
    if hgnc_ens is None:
        hgnc_ens = '{}/data/hgncid_name_ens.txt'.format(curdir)

    # Get list of allowed genes
    cmd1 = [
        '{}/scripts/get_allowlist.py'.format(curdir),
        '--gtf', options.gencode, 
        '-he', hgnc_ens,
        '-o', options.output_prefix, 
        '--mito', options.mito_name_hg38
        ]
    print(" ".join(cmd1), file=sys.stderr)
    subprocess.call(cmd1)
    
    if options.human_annotation is not None:
        # Filter the corresponding human annotation before we delete
        # necessary tmp files
        cmd1h = [
            '{}/scripts/filter_human_annotation.py'.format(curdir),
            '--gtf', options.human_annotation,
            '--allowed_genes', '{}.genes'.format(options.output_prefix),
            '--allowed_tx', '{}.tx'.format(options.output_prefix),
            '--tx2gene', '{}.tx2gene'.format(options.output_prefix),
            '--mito', '{}.mito'.format(options.output_prefix),
            '-o', '{}.human.unsorted.gtf'.format(options.output_prefix)
        ]
        print(" ".join(cmd1h), file=sys.stderr)
        subprocess.call(cmd1h)

        cmd2h = [
            '{}/scripts/sort_annotation.py'.format(curdir),
            '-i', '{}.human.unsorted.gtf'.format(options.output_prefix),
            '-o', '{}.human.gtf'.format(options.output_prefix)
        ]
        print(" ".join(cmd2h), file=sys.stderr)
        subprocess.call(cmd2h)

        os.unlink('{}.human.unsorted.gtf'.format(options.output_prefix))


    # Filter the CAT annotation and remove mitochondrial genes
    cmd2 = [
        '{}/scripts/cleanup_cat_annotation.py'.format(curdir), 
        '--gff3', options.cat_gff3, 
        '-a', options.output_prefix, 
        '-o', '{}'.format(options.output_prefix)
    ]
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
    cmd3 = [
        '{}/scripts/mito_liftoff.py'.format(curdir), 
        '--cat_fasta', options.cat_fasta, 
        '--hg38_fasta', options.hg38_fasta,
        '--cat_mito', options.cat_mito, 
        '--allowlist_base', options.output_prefix, 
        '-o', options.output_prefix
    ]
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

    if options.scavenge_ensembl_annotation:
        
        if os.path.isfile('{}/data/{}_ids.txt'.format(curdir, options.species_id)):
            # It works
            pass
        else:
            # Check that the user has given a valid name key for biomaRt
            p = subprocess.Popen([
                '{}/scripts/get_ensembl_names.R'.format(curdir)
            ], stdout=subprocess.PIPE)
            out, err = p.communicate()
            found_id = False
            for line in out.decode().split('\n'):
                dat = line.split('\t')
                if dat[0] == options.species_id:
                    found_id = True
                    break
            if not found_id:
                print("ERROR: species id {} not found in Ensembl.".\
                    format(options.species_id), file=sys.stderr)
                print("To see a list of eligible IDs, please run \
{}/scripts/get_ensembl_names.R.".format(curdir), file=sys.stderr)
                print("NOTE: it's also possible that a curl query timed \
out, since that seems to happen a lot with Ensembl web servers. Check \
output for error messages and run again if that's the case.", file=sys.stderr)
                exit(1)

        # Attempt to rescue by lifting from an old Ensembl annotation
        cmd4 = [
            '{}/scripts/get_homologous_genelist.R'.format(curdir), 
            '{}.rm.genes.notx'.format(options.output_prefix),
            options.species_id
        ]
        f = open('{}.homology'.format(options.output_prefix), 'w')
        print(" ".join(cmd4), file=sys.stderr)
        p = subprocess.Popen(cmd4, stdout=f)
        out, err = p.communicate()
        f.close()
        
        # Check that the user spelled the name right
        n_hom_lines = 0
        f = open('{}.homology'.format(options.output_prefix), 'r')
        for line in f:
            n_hom_lines += 1
        f.close()
        if n_hom_lines == 0:
            print("WARNING: no gene mappings found. There are no genes in \
Ensembl homologous to the ones that were eliminated.".\
                format(options.species_id), file=sys.stderr)
        else:    
            if options.chain is None:
                # Same genome assembly
                
                # Get a file mapping chromosome names between assemblies
                if not os.path.isfile('{}.fai'.format(options.ensembl_fasta)):
                    print("Indexing {}...".format(options.ensembl_fasta), file=sys.stderr)
                    subprocess.call(['samtools', 'faidx', options.ensembl_fasta])
                if not os.path.isfile('{}.fai'.format(options.cat_fasta)):
                    print("Indexing {}...".format(options.cat_fasta), file=sys.stderr)
                    subprocess.call(['samtools', 'faidx', options.cat_fasta])
                
                outf = open('{}.ensmap'.format(options.output_prefix), 'w')

                cmd5 = [
                    '{}/scripts/map_between_assemblies.py'.format(curdir),
                    '{}.fai'.format(options.ensembl_fasta),
                    '{}.fai'.format(options.cat_fasta)
                ]

                print(" ".join(cmd5), file=sys.stderr)
                p = subprocess.Popen(cmd5, stdout=outf)
                out, err = p.communicate()
                outf.close()

                cmd6 = [
                    '{}/scripts/scavenge_old_ensembl_annotation.py'.format(curdir),
                    '-G', '{}.homology'.format(options.output_prefix),
                    '-g', options.old_ens_gtf,
                    '--hgnc_ens', hgnc_ens,
                    '--idmap', '{}.ensmap'.format(options.output_prefix),
                    '-o', options.output_prefix
                ]

                print(" ".join(cmd6), file=sys.stderr)
                subprocess.call(cmd6)
                
                os.unlink('{}.ensmap'.format(options.output_prefix))

            else:
                # Different genome assembly
                
                # Get a file mapping chromosome names between assemblies
                if not os.path.isfile('{}.fai'.format(options.ucsc_fasta)):
                    print("Indexing {}...".format(options.ucsc_fasta), file=sys.stderr)
                    subprocess.call(['samtools', 'faidx', options.ucsc_fasta])
                if not os.path.isfile('{}.fai'.format(options.cat_fasta)):
                    print("Indexing {}...".format(options.cat_fasta), file=sys.stderr)
                    subprocess.call(['samtools', 'faidx', options.cat_fasta])
                
                outf = open('{}.ucscmap'.format(options.output_prefix), 'w')

                cmd5 = [
                    '{}/scripts/map_between_assemblies.py'.format(curdir),
                    '{}.fai'.format(options.ucsc_fasta),
                    '{}.fai'.format(options.cat_fasta)
                ]

                print(" ".join(cmd5), file=sys.stderr)
                p = subprocess.Popen(cmd5, stdout=outf)
                out, err = p.communicate()
                outf.close()

                cmd6 = [
                    '{}/scripts/lift_old_ensembl_annotation.py'.format(curdir), 
                    '-G', '{}.homology'.format(options.output_prefix),
                    '-g', options.old_ens_gtf, 
                    '-c', options.chain, 
                    '--hgnc_ens', hgnc_ens, 
                    '--idmap', '{}.ucscmap'.format(options.output_prefix),
                    '-o', options.output_prefix
                ]
                if options.gtfToGenePred_path is not None:
                    cmd6.append("--gtfToGenePred_path")
                    cmd6.append(options.gtfToGenePred_path)
                if options.liftOver_path is not None:
                    cmd6.append("--liftOver_path")
                    cmd6.append(options.liftOver_path)
                if options.genePredToGtf_path is not None:
                    cmd6.append("--genePredToGtf_path")
                    cmd6.append(options.genePredToGtf_path)

                print(" ".join(cmd6), file=sys.stderr)
                subprocess.call(cmd6)
                
                os.unlink('{}.ucscmap'.format(options.output_prefix))
            
            os.unlink('{}.homology'.format(options.output_prefix))

            # Alter removed gene list
            gidrescue = set([])
            f = open('{}.rescued.gid'.format(options.output_prefix), 'r')
            for line in f:
                line = line.rstrip()
                if line != "":
                    gidrescue.add(line)
            f.close()
            os.unlink('{}.rescued.gid'.format(options.output_prefix))

            gidrm = []
            f = open('{}.rm.genes'.format(options.output_prefix), 'r')
            for line in f:
                line = line.rstrip()
                dat = line.split('\t')
                if len(dat) >= 3 and dat[2] != "":
                    if dat[2] not in gidrescue:
                        gidrm.append(line)
                else:
                    gidrm.append(line)
            f.close()
            f = open('{}.rm.genes'.format(options.output_prefix), 'w')
            for line in gidrm:
                print(line, file=f)
            f.close()

    # Clean up
    if os.path.isfile('{}.tx'.format(options.output_prefix)):
        os.unlink('{}.tx'.format(options.output_prefix))
    if os.path.isfile('{}.tx2gene'.format(options.output_prefix)):
        os.unlink('{}.tx2gene'.format(options.output_prefix))
    if os.path.isfile('{}.genes'.format(options.output_prefix)):
        os.unlink('{}.genes'.format(options.output_prefix))
    if os.path.isfile('{}.rm.tx'.format(options.output_prefix)):
        os.unlink('{}.rm.tx'.format(options.output_prefix))
    if os.path.isfile('{}.rm.genes.notx'.format(options.output_prefix)):
        os.unlink('{}.rm.genes.notx'.format(options.output_prefix))

    # Now put everything together into one annotation.
    cmd6 = [
        "{}/scripts/sort_annotation.py".format(curdir),
        "-o", "{}.gtf".format(options.output_prefix),
        "-i", "{}.main.gtf".format(options.output_prefix),
        "{}.mito.liftoff.gtf".format(options.output_prefix)
    ]

    if options.scavenge_ensembl_annotation:
        if os.path.isfile('{}.rescued.gtf'.format(options.output_prefix)):
            cmd6.append("{}.rescued.gtf".format(options.output_prefix))

    print(" ".join(cmd6), file=sys.stderr)
    subprocess.call(cmd6)

    # Clean up
    os.unlink("{}.main.gtf".format(options.output_prefix))
    os.unlink("{}.mito.liftoff.gtf".format(options.output_prefix))
    if os.path.isfile('{}.rescued.gtf'.format(options.output_prefix)):
        os.unlink("{}.rescued.gtf".format(options.output_prefix))


if __name__ == '__main__':
    sys.exit(main(sys.argv))
