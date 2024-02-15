# litterbox
Clean up CAT (Comparative Annotation Toolkit) annotations

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

#### Clean up the CAT GTF - autosomal annotations

1. Run `gff32gtf.py` on the CAT GFF3 file. This changes the sort order and the names of things in the input GFF3 to make `cellranger mkref` happy. It then outputs GTF, which is used by popular indexing/alignment programs.
2. Run `get_allowlist.py` on the most recent [GENCODE GTF](https://www.gencodegenes.org/human/) and a tab separated 3-column file mapping `HGNC ID  Approved Symbol  Ensembl Gene ID` from the [HUGO Gene Nomenclature Committee](https://www.genenames.org/download/custom/). This outputs a filtered list of ENSEMBL genes and approved symbol for genes that pass the [10X Genomics filtering criteria](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps).
3. Run `get_allowlist.py` again with same data and the `-m` option to output a similar list, for mitochondrial genes only.
4. Run `cleanup_cat_annotation.py` with the GTF output by `gff32gtf.py`, the allow list created by `get_allowlist.py`, and the list of mitochondrial genes output by `get_allowlist.py`. The result is a filtered GTF without mitochondrial gene annotations.

 #### Make new, clean mitochondrial gene annotations

 1. Pull all mitochondrial annotations from the latest gencode, using the list of mitochondrial genes you created in step 3 above. If the file is called `mito.txt` and the annotation is called `gencode.gtf.gz`, then run `zcat gencode.gtf.gz | grep -f -w mito.txt > gencode.mito.gtf`
 2. Get the mitochondrial genome from hg38 (assuming this is still the current genome version that Gencode is mapped to)
 3. Get the mitochondrial genome for your organism. If your FASTA file has been indexed with samtools faidx, this can be done like `samtools faidx genome.fa chrM > chrM.fa`
 4. Run [liftOff](https://github.com/agshumate/Liftoff) to lift the mitochondrial GTF from the hg38 mitochondrial sequence onto your organism
 5. Add these mitochondrial annotations to your annotation from the section above using `cat`

#### Make a new version of the human annotation with the names you used in the CAT annotation
 1. Run `alter_gene_names_human.py` using your allow list from the first section and your human GTF. Outputs a new human GTF

#### Index everything and run your alignments

Giddy up


