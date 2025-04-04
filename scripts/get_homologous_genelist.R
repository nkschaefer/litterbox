#! /usr/bin/env Rscript

# If some genes were lost from a CAT annotation, give the list of
# genes here (as Human Ensembl IDS), plus the species name
# as the lowercase first letter of species + lowercase genus,
# no spaces.

# Finds the homologous genes in the other species.


args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
    write("Please provide *.rm.genes.notx and homology-mapping file", stderr())
    q()
}

gfile <- read.table(args[1], sep="\t")
# ENSG is third
gfile <- gfile[which(gfile$V3 != ""),]

genelist <- unique(gfile$V3)

# Check to see if we already downloaded IDs.
human2other <- read.table(args[2], sep='\t')
human2other <- human2other[which(human2other$V1 %in% genelist),]
write.table(human2other, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

