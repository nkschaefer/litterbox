# litterbox
Clean up CAT (Comparative Annotation Toolkit) annotations

## Requirements

* [samtools](https://github.com/samtools/samtools)
* gtfToGenePred, genePredToGtf, and liftOver from [UCSC Kent utils](https://hgdownload.soe.ucsc.edu/admin/exe/)
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
* [pandas](https://pandas.pydata.org/)
* [natsort](https://pypi.org/project/natsort/)
* [liftoff](https://github.com/agshumate/Liftoff)

## Installation

The easiest thing to do is to use [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to install all the requirements using the included yml file:

`conda env create --file litterbox.yml`

## Why

[CAT](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit) makes good annotations for cross-species analysis. But its GFF3 files can pose problems for cross-species DE analysis.

### What's wrong

* Mitochondrial annotations tend to be messed up
* Mitochondrial genes can appear on the autosomes, especially where there are nuclear mitochondrial insertions (NUMTs)
* Source annotations can include features that steal reads from other genes and confound RNA-seq quantification, including [read-through transcripts](https://www.ensembl.info/2019/02/11/annotating-readthrough-transcription-in-ensembl/) that are [filtered out](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps) by default in the 10X mkref pipeline.
* De novo annotations can be wrong and cause the same problem
* GFF3s are not necessarily sorted and formatted in the way expected by indexing programs (like `STAR genomeGenerate` and `cellranger mkref`)
* When using tools like [scanpy](https://scanpy.readthedocs.io/en/stable/) or [Seurat](https://satijalab.org/seurat/), genes are usually keyed to their human-readable name, not their unique IDs. If you are using annotations from multiple species, [gene synonyms](https://www.genenames.org/tools/multi-symbol-checker/) might mean that genes are dropped when they don't have the same name for one or more species. 

### How to fix it

#### Get data

`litterbox` needs the following things:

* The CAT annotation you want to clean up (GFF3 format)
* A recent [Gencode](https://www.gencodegenes.org/human/) (Human) annotation (GTF format)
* The name of the mitochondrial sequence in the genome assembly the CAT annotation corresponds to (it's probably chrM; you can check by indexing the FASTA using `samtools faidx` and looking at the sequence names in the first column of the resulting `fai` file
* The genome (FASTA format) the HAL annotation corresponds to. If you have a set of aligned genomes in [HAL](https://github.com/ComparativeGenomicsToolkit/hal) format, you can extract the genome you want with the HAL toolkit command `hal2fasta` followed by the name of the genome in the HAL file.
* The genome (FASTA format) that Gencode corresponds to: this is [hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
* a tab-separated file, mapping (in three columns) `HGNC ID` -> `Approved symbol` -> `Ensembl Gene ID`. You can generate an up-to-date one [here](https://www.genenames.org/download/custom/) or use the one included in the `data` directory (default).

Recommended other things, to rescue some of the genes that will otherwise get dropped:

`litterbox` can also try to pull in records for dropped genes from an [Ensembl](https://www.ensembl.org) annotation. This is a last resort: if all transcripts that CAT used to place a gene get filtered out by their tags/biotypes in the latest Gencode annotation, all we can do is try to pull some of those genes back from a different annotation. There are two possibilities:

#### If Ensembl annotated the most recent assembly

If there is an Ensembl annotation of the same organism on the same genome assembly version (you can check dates against those in the UCSC genome browser, or if still uncertain, compare sequence lengths in the .fai files), then you will need to download the GTF of the Ensembl annotation and the genome assembly the annotation used, in FASTA format. 

You will also need to know the name of the organism of the CAT assembly as it exists in the Ensembl database. This should be lowercase first letter of species name + lowercase genus name, with no spaces. Example: hsapiens

#### If Ensembl only annotated an older assembly

If the most recent Ensembl annotation for the organism is on a different assembly version than the one that CAT annotated, then you will need to download that assembly (in GTF format), the UCSC [chain file](https://genome.ucsc.edu/goldenPath/help/chain.html) for lifting from that assembly version to the one used for CAT ([UCSC downloads](https://hgdownload2.soe.ucsc.edu/downloads.html)), and a FASTA file of the UCSC version of the same assembly as was used for CAT. 

You will also need to know the name of the organism of the CAT assembly as it exists in the Ensembl database. This should be lowercase first letter of species name + lowercase genus name, with no spaces. Example: hsapiens

### Run the program

Everything is set up so you can just run one program (in the main directory): `filter_annotation.py`. It has a lot of arguments that will show up if you run it with `--help` or `-h.` You must provide an output file name prefix, and the final annotation will be in `[output_prefix].gtf`. 

You can optionally filter a human Gencode annotation (recommended: the version of Gencode that was used to create the CAT annotation) to contain the same genes, and the same gene names, as the final CAT annotation. This is a good idea to make sure that gene names match up and that genes that exist in both annotations, but left included in only one, don't drive up/down expression of other genes as a result of being left in or out. If you choose this option, the resulting annotation will be `[output_prefix].human.gtf`.

### What it's doing

1. Goes through the latest Gencode GTF and makes a list of allowed gene IDs and allowed transcript IDs, based on values of `gene_type`, `transcript_type`, and `tag` fields. Allowed type fields were taken from the 10X Genomics `cellranger mkref` instructions: `protein_coding`, `lncRNA`, `IG_C_gene`, `IG_D_gene`, `IG_LV_gene`, `IG_V_gene`, `IG_V_pseudogene`, `IG_J_pseudogene`, `IG_C_pseudogene`, `TR_C_gene`, `TR_D_gene`, `TR_J_gene`, `TR_V_gene`, `TR_V_pseudogene`, and `TR_J_pseudogene`. Disallowed tag fields were also based on this guide, with the addition of several more: `readthrough_gene`, `readthrough_transcript`, `PAR`, `fragmented_locus`, and `low_sequence_quality`. 
2. Goes through the CAT GFF3, converts to GTF, and removes all genes and transcripts that were lifted from the human gencode annotation and did not pass the filters in step 1. If a transcript was a prediction not based on homology with human (i.e. the Gencode annotation), it's removed if `transcript_class` is `poor_alignment`, or if `transcript_class` is `possible_paralog` and either `reference_support` is not true, or both `rna_support` and `pacbio_isoform_supported` are not true. If all transcripts of a gene are removed, then the gene is also removed.
3. Renames transcripts that pass filters to their preferred symbol/common name, according to the latest HGNC data (included in `data` directory)
4. Removes all mitochondrial genes (as determined from the recent Gencode annotation) from the CAT annotation, to avoid cases where [NUMTs](https://en.wikipedia.org/wiki/Nuclear_mitochondrial_DNA_segment) are annotated as mitochondrial genes instead of the actual mitochondrial genes
5. Uses [liftoff](https://github.com/agshumate/Liftoff) and the mitochondrial sequences from hg38 and the CAT-annotated assembly to re-annotate the mitochondrial genome in the CAT-annotated assembly
6. Optionally attempts to rescue some of the eliminated genes by taking their annotations from an Ensembl gene annotation. First, takes the list of genes that were eliminated because of all their transcripts being eliminated (and not directly filtered out themselves). Then converts these human Ensembl gene IDs to IDs for the species of interest, using [biomaRt](https://www.ensembl.org/info/data/biomart/index.html).
    1. If the Ensembl gene annotation you downloaded is on the same assembly coordinates, pulls the relevant genes out of the annotation and converts sequence names to those used in CAT by mapping chromosome names to each other across assemblies (matching them up based on their sequence lengths in `fai` files for both genomes).
    2. If the Ensembl gene annotation you downloaded is on different/older assembly coordinates, pulls the relevant genes out of the annotation and then lifts this annotation to the correct assembly, using the UCSC Kent tools `gtfToGenePred`, `liftOver`, and `genePredToGtf`. Then converts sequence names to those used in CAT by mapping chromosome names to each other across assemblies (matching them up based on their sequence lengths in `fai` files for both genomes). For a gene to survive make it into the final annotation, it must have a `gene` feature and at least one `transcript` and at least one `exon` feature on the new genome coordinates.
 7. Combines the new main GTF, the mitochondrial GTF, and (optionally) the "rescued" GTF and sorts these in a way designed to make `cellranger-arc mkref` happy. The sort order matters to this program, as documented [here](https://github.com/10XGenomics/cellranger/issues/133), and the code to sort properly was taken from [chbk's suggestion](https://github.com/10XGenomics/cellranger/issues/133#issuecomment-989119662) in that issue (it's in the repository in `scripts/sort_annotation.py`).
 8. Optionally also creates a new human annotation with the same filters: if you provide a human annotation (such as Gencode of the same version used to make the CAT annotation), filters the genes the same way as the CAT annotation, and sets the `gene_name` field to the HGNC approved symbol for each gene, provided one exists. If you are cleaning up a set of CAT annotations for multiple species all based on the same Gencode, you only need to do this once to have a corresponding human annotation. The human annotation is adjusted according to the HGNC IDs and the information in the latest Gencode release, but does not depend on the provided CAT annotation.

### How about that

You got yourself a clean annotation


