#! /usr/bin/env Rscript

# If some genes were lost from a CAT annotation, give the list of
# genes here (as Human Ensembl IDS), plus the species name
# as the lowercase first letter of species + lowercase genus,
# no spaces.

# Finds the homologous genes in the other species.


args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2){
    write("Please provide *.rm.genes.notx and species identifier", stderr())
    write("Example identifiers: ptroglodytes, mmulatta", stderr())
    q()
}

gfile <- read.table(args[1], sep="\t")
# ENSG is third
gfile <- gfile[which(gfile$V3 != ""),]

genelist <- unique(gfile$V3)

# Check to see if we already downloaded IDs.
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
dir.root <- dirname(dirname(script.name))
fn <- paste(dir.root, '/data/', args[2], '_ids.txt', sep='')
if (file.exists(fn)){
    human2other <- read.table(fn)
    human2other <- human2other[which(human2other$V1 %in% genelist),]
    write.table(human2other, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    q()
}

#library(httr)
#setconfig(config(sslverifypeer = 0L))

suppressPackageStartupMessages(library(biomaRt))

human <- useEnsembl('genes', 'hsapiens_gene_ensembl')

field_id <- paste(args[2], '_homolog_ensembl_gene', sep='')

# Get a list of homologous gene IDs in that organism's genome
res <- getBM(c("ensembl_gene_id", field_id), filters=c("ensembl_gene_id"), values=genelist, mart=human)
res <- res[which(res[[field_id]] != ""),]

write.table(res, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

