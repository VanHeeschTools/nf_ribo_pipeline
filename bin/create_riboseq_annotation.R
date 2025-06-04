#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(RiboseQC)
  library(rtracklayer)
})

# Set global variables ----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

twobit_file <- args[1]
gtf <- args[2]
package_install_loc <- args[3]
annot_name <- "test_annotation_package"

# Current problem: R tries to use a tmp directory that is not writable (container dir)
# I tried this to make sure R knows what it can use, but it did not work unfortunately.
# The biggest problem is that this script works outside nextflow, so the problem could also
# lie there.
Sys.setenv(TMPDIR=package_install_loc)
configure.vars=paste0("TMPDIR=", package_install_loc)

paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Prepare annotation files ------------------------------------------------
message("Preparing annotation ...")
prepare_annotation_files(annotation_directory = package_install_loc,
                         twobit_file = twobit_file,
                         gtf_file = gtf,
                         annotation_name = annot_name,
                         forge_BSgenome = TRUE)

BSgenome_dir <- grep("BSgenome", x = list.dirs(package_install_loc,
                                               recursive = FALSE),
                     value = TRUE)

message("Installing annotation ...")
install.packages(BSgenome_dir,
                 character.only = TRUE,
                 repos = NULL,
                 type = "source")
