#!/usr/bin/env Rscript

### Run this script to generate a mapping of genes between species



### user defined variables

# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE)
args <- strsplit(cli, "=", fixed = TRUE)
args

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}



### packages
library(dplyr)
library(tidyr)
library(GenomicFeatures)
library(txdbmaker)
library(rtracklayer)
library(BSgenome)
library(gprofiler2)
library(biomaRt)

# Check and install the species-specific BSgenome packages
# The variable is a string, so we must load it using character.only = TRUE
library(BiocManager)
if (!require(biomart_genome, character.only = TRUE)) BiocManager::install(biomart_genome)
if (!require(human_biomart_genome, character.only = TRUE)) BiocManager::install(human_biomart_genome)





### functions

#' @title Generate Transcriptome FASTA (Exon-based)
#' @description This function retrieves exon coordinates for a specific organism from BioMart, 
#' creates a TranscriptDb object, collapses overlapping exons per gene (creating a flattened 
#' gene model), extracts the DNA sequence from the BSgenome object, maps Ensembl IDs to 
#' external gene symbols, and writes the resulting sequences to a FASTA file.
#'
#' @param organism Character. The organism name (used for filtering/printing).
#' @param chrom_number Integer. The number of autosomes to include (e.g., 22 for humans). 
#' Function also automatically attempts to include X and Y if present.
#' @param dataset Character. The Ensembl BioMart dataset name (e.g., "hsapiens_gene_ensembl").
#' @param genome Character. The name of the BSgenome package/object to extract sequences from.
#' @param out_dir Character. The directory path where the output FASTA file will be saved.
#' @param out_name Character. The prefix for the output filename (suffix will be "_transcriptome.fa").
#'
#' @return A DNAStringSet object containing the extracted sequences, invisibly. 
#' Side effect: Writes a file named `{out_name}_transcriptome.fa` to `out_dir`.
make_transcriptome_fasta <- function(organism, chrom_number, dataset, genome, out_dir, out_name) {
  
  # create transcript database
  txdb <- txdbmaker::makeTxDbFromBiomart(dataset=dataset)

  # subset to main chromosomes
  chr_keep <- as.character(seq(1, chrom_number))

  if ("X" %in% seqlevels(txdb)) {
    chr_keep <- c(chr_keep, "X")
  }
  if ("Y" %in% seqlevels(txdb)) {
    chr_keep <- c(chr_keep, "Y")
  }

  seqlevels(txdb) <- chr_keep

  # convert chromosome names to UCSC format
  newSeqNames <- paste('chr', seqlevels(txdb), sep = '')
  names(newSeqNames) <- seqlevels(txdb)
  txdb <- renameSeqlevels(txdb, newSeqNames)

  # generate GRanges list object, listing exons per gene
  tmp <- exonsBy(txdb, by="gene")

  # combine overlapping exons for individual genes
  tmp1 <- reduce(tmp)

  # get DNA sequences
  fasta <- extractTranscriptSeqs(eval(as.name(genome)), tmp1)

  # convert gene IDs to gene symbols
  #ids2names <- gconvert(query = fasta@ranges@NAMES, organism = organism, 
  #                      target="ENSG", numeric_ns="ENTREZGENE_ACC", 
  #                      mthreshold = Inf, filter_na = TRUE)
  #ids2names <- ids2names[,c("input","name")]

  ensembl <- useEnsembl(biomart="ensembl", dataset=dataset)
  ids2names <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = fasta@ranges@NAMES, mart = ensembl)
  colnames(ids2names) <- c("input", "name")

  ids2names <- ids2names[ids2names$name!="None",]
  ids2names <- ids2names[!duplicated(ids2names$input),]

  # find gene symbols of gene IDs in sequence object
  fasta_names <- as.data.frame(fasta@ranges@NAMES)
  colnames(fasta_names) <- "gene_id"
  fasta_id_names <- fasta_names %>% left_join(ids2names[,c("input","name")], 
                                              by=join_by(gene_id==input))

  # check that order of gene ID to gene symbol file matches that of gene IDs in sequence object
  print(paste("Name order matches: ", identical(fasta_id_names$gene_id, fasta@ranges@NAMES), sep=""))

  # replace gene IDs in sequence object with gene symbols
  fasta@ranges@NAMES <- fasta_id_names$name

  # drop genes without names
  fasta <- fasta[!is.na(fasta@ranges@NAMES),]

  entries_keep <- list()
  for (i in seq(1, length(fasta@ranges@NAMES))) {
    if (nchar(fasta@ranges@NAMES[i]) > 0) {
      entries_keep <- append(entries_keep, i)
    }
  }
  entries_keep <- unlist(entries_keep)
  fasta <- fasta[entries_keep,]

  # save as fasta file
  writeXStringSet(fasta, filepath=paste(out_dir, "/", out_name, "_transcriptome.fa", sep=""),
                  format="fasta")
  return(fasta)             
}



#' @title Generate Coding Sequence (CDS) FASTA
#' @description This function retrieves Coding Sequence (CDS) coordinates for a specific organism 
#' from BioMart, creates a TranscriptDb object, collapses overlapping CDS regions per gene, 
#' extracts the DNA sequence from the BSgenome object, maps Ensembl IDs to external gene symbols, 
#' and writes the resulting sequences to a FASTA file.
#'
#' @param organism Character. The organism name (used for filtering/printing).
#' @param chrom_number Integer. The number of autosomes to include.
#' @param dataset Character. The Ensembl BioMart dataset name (e.g., "hsapiens_gene_ensembl").
#' @param genome Character. The name of the BSgenome package/object to extract sequences from.
#' @param out_dir Character. The directory path where the output FASTA file will be saved.
#' @param out_name Character. The prefix for the output filename (suffix will be "_cds.fa").
#'
#' @return A DNAStringSet object containing the extracted CDS sequences, invisibly. 
#' Side effect: Writes a file named `{out_name}_cds.fa` to `out_dir`.
make_cds_fasta <- function(organism, chrom_number, dataset, genome, out_dir, out_name) {
  
  # create transcript database
  txdb <- txdbmaker::makeTxDbFromBiomart(dataset=dataset)

  # subset to main chromosomes
  chr_keep <- as.character(seq(1, chrom_number))

  if ("X" %in% seqlevels(txdb)) {
    chr_keep <- c(chr_keep, "X")
  }
  if ("Y" %in% seqlevels(txdb)) {
    chr_keep <- c(chr_keep, "Y")
  }

  seqlevels(txdb) <- chr_keep

  # convert chromosome names to UCSC format
  newSeqNames <- paste('chr', seqlevels(txdb), sep = '')
  names(newSeqNames) <- seqlevels(txdb)
  txdb <- renameSeqlevels(txdb, newSeqNames)

  print(seqlevels(txdb))

  # generate GRanges list object, listing coding sequences per gene
  tmp <- cdsBy(txdb, by="gene")

  # combine overlapping exons for individual genes
  tmp1 <- reduce(tmp)

  # get DNA sequences
  fasta <- extractTranscriptSeqs(eval(as.name(genome)), tmp1)

  # convert gene IDs to gene symbols
  #ids2names <- gconvert(query = fasta@ranges@NAMES, organism = organism, 
  #                      target="ENSG", numeric_ns="ENTREZGENE_ACC", 
  #                      mthreshold = Inf, filter_na = TRUE)
  #ids2names <- ids2names[,c("input","name")]

  ensembl <- useEnsembl(biomart="ensembl", dataset=dataset)
  ids2names <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters = 'ensembl_gene_id', values = fasta@ranges@NAMES, mart = ensembl)
  colnames(ids2names) <- c("input", "name")

  ids2names <- ids2names[ids2names$name!="None",]
  ids2names <- ids2names[!duplicated(ids2names$input),]

  # find gene symbols of gene IDs in sequence object
  fasta_names <- as.data.frame(fasta@ranges@NAMES)
  colnames(fasta_names) <- "gene_id"
  fasta_id_names <- fasta_names %>% left_join(ids2names[,c("input","name")], 
                                              by=join_by(gene_id==input))

  # check that order of gene ID to gene symbol file matches that of gene IDs in sequence object
  print(paste("Name order matches: ", identical(fasta_id_names$gene_id, fasta@ranges@NAMES), sep=""))

  # replace gene IDs in sequence object with gene symbols
  fasta@ranges@NAMES <- fasta_id_names$name

  # drop genes without names
  fasta <- fasta[!is.na(fasta@ranges@NAMES),]

  entries_keep <- list()
  for (i in seq(1, length(fasta@ranges@NAMES))) {
    if (nchar(fasta@ranges@NAMES[i]) > 0) {
      entries_keep <- append(entries_keep, i)
    }
  }
  entries_keep <- unlist(entries_keep)
  fasta <- fasta[entries_keep,]

  # save as fasta file
  writeXStringSet(fasta, filepath=paste(out_dir, "/", out_name, "_cds.fa", sep=""),
                  format="fasta")
  return(fasta)             
}



### Analysis
print(paste("Generating coding sequence fasta file for ", species, sep="")) 

if (!file.exists(paste(out_dir, "/", species, "_cds.fa", sep=""))) {
  make_cds_fasta(organism=species,
                chrom_number=chrom_number,
                dataset=biomart_dataset, 
                genome=biomart_genome,
                out_dir=out_dir, 
                out_name=species)
}

print("Generating coding sequence fasta file for human") 

if (!file.exists(paste(out_dir, "/human_cds.fa", sep=""))) {
  make_cds_fasta(organism="human",
                chrom_number=22,
                dataset=human_biomart_dataset, 
                genome=human_biomart_genome,
                out_dir=out_dir, 
                out_name="human")
}