# litterbox
Clean up CAT (Comparative Annotation Toolkit) annotations

## Why

[CAT](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit) is a powerful tool that can make high-quality annotations for cross-species analysis. But its GFF3 files can pose a few problems for cross-species DE analysis.

## What's wrong

* Mitochondrial annotations tend to be messed up
* Mitochondrial genes can appear on the autosomes, especially where there are nuclear mitochondrial insertions (NUMTs)
* Source annotations can include features that steal reads from other genes and confound RNA-seq quantification, including [read-through transcripts](https://www.ensembl.info/2019/02/11/annotating-readthrough-transcription-in-ensembl/) that are [filtered out](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps) by default in the 10X mkref pipeline.
* De novo annotations can be wrong and cause the same problem
* GFF3s are not necessarily sorted and formatted in the way expected by indexing programs (like S`TAR genomeGenerate` and `cellranger mkref`)

## How to fix it

* `gff32gtf.py` changes the sort order and the names of things in the input GFF3 to make `cellranger mkref` happy. It then outputs GTF, which is used by popular indexing/alignment programs.
* 

