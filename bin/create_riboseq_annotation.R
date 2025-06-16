#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(RiboseQC)
  library(BSgenome)
  library(Biostrings)
})

# Set global variables ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

twobit_file <- args[1]
gtf <- args[2]
genome <- args[3]

annot_name <- "custom"
seed_file <- file.path(getwd(), paste0(annot_name, ".seed"))

Sys.setenv(TMPDIR=getwd())
configure.vars=paste0("TMPDIR=", getwd())
paths <- c(getwd(), .libPaths())
.libPaths(paths)

# Prepare annotation files ------------------------------------------------
message("Preparing annotation ...")
prepare_annotation_files(annotation_directory = getwd(),
                         twobit_file = twobit_file,
                         gtf_file = gtf,
                         genome_seq = genome,
                         annotation_name = annot_name,
                         forge_BSgenome = F)


message("Creating custom BSgenome package using: ")
message(seed_file)

seed_lines <- c(
  "Package: BSgenome.Hsapiens.nextflow",
  "Title: Custom BSgenome Package",
  "Description: Built from multi-FASTA",
  "Version: 1.0.0",
  "organism: Homo_sapiens",
  "common_name: Human",
  "provider: Custom",
  "provider_version: Custom_v1",
  "BSgenomeObjname: Hsapiens",
  "genome: Custom_nextflow",
  paste0("seqs_srcdir: ", dirname(twobit_file),"/"),
  paste0("seqfile_name: ", basename(twobit_file)),
  "seqnames: c(1:22, 'X', 'Y')",
  "circ_seqs: character(0)",
  "organism_biocview: Homo_sapiens"
)

writeLines(seed_lines, con = seed_file)


forgeBSgenomeDataPkg(
  seed_file,        # Path to the seed file
  destdir = savedir  # Output location
)

message("Forged BSgenome package, installing package")

BSgenome_dir <- grep("BSgenome", x = list.dirs(getwd(),
                                                recursive = FALSE),
                      value = TRUE)

install.packages(getwd(),
                  character.only = TRUE,
                  repos = NULL,
                  type = "source")


message("Finished installing custom BSgenome package")