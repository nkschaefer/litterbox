#! /usr/bin/env Rscript

library(biomaRt)

human <- useEnsembl('genes', 'hsapiens_gene_ensembl')

attr <- listAttributes(human)
attr <- attr[which(attr$page=="homologs"),]
attr <- attr[grep("_homolog_ensembl_gene$", attr$name),c(1,2)]
attr$description <- gsub(" gene stable ID", "", attr$description)
attr$name <- gsub("_homolog_ensembl_gene", "", attr$name)
attr <- attr[order(attr$name),]

write.table(attr, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

