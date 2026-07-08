#!/usr/bin/env Rscript

#' create_riboseq_annotation.R
#'
#' Prepares RiboseQC annotation files (BSgenome package + TxDb) from a
#' reference genome (2bit and fa) and GTF file.
#'
#' Usage:
#'   create_riboseq_annotation.R <twobit_file> <gtf> <package_install_loc> \
#'                                <annot_name> <savedir> <genome>

report_error_and_exit <- function(stage, err) {
    message("")
    message("ERROR during stage: ", stage)
    message(conditionMessage(err))
    message("")
    message("Either check if the input files are all correct, or run the script outside the pipeline and provide input parameters for:")
    message("package_install_loc")
    message("orfquant_annotation")
    quit(status = 1, save = "no")
}

# Parse and validate arguments 
tryCatch({
    args <- commandArgs(trailingOnly = TRUE)

    if (length(args) < 6) {
        stop(sprintf(
            "Expected 6 arguments (twobit_file, gtf, package_install_loc, annot_name, savedir, genome), got %d: %s",
            length(args), paste(args, collapse = ", ")
        ))
    }

    twobit_file         <- args[1]
    gtf                 <- args[2]
    package_install_loc <- args[3]
    annot_name          <- args[4]
    savedir             <- args[5]
    genome              <- args[6]

    message("Arguments received:")
    message("  twobit_file          = ", twobit_file)
    message("  gtf                  = ", gtf)
    message("  package_install_loc  = ", package_install_loc)
    message("  annot_name           = ", annot_name)
    message("  savedir              = ", savedir)
    message("  genome               = ", genome)

    # Basic input file existence checks - fail fast with a clear message
    # rather than letting a downstream package throw a cryptic error.
    if (!file.exists(twobit_file)) {
        stop("2bit file not found: ", twobit_file)
    }
    if (!file.exists(gtf)) {
        stop("GTF file not found: ", gtf)
    }
}, error = function(e) report_error_and_exit("argument parsing / validation", e))

# Set up TMPDIR and library paths 
tryCatch({
    Sys.setenv(TMPDIR = package_install_loc)
    configure.vars <- paste0("TMPDIR=", package_install_loc)

    paths <- c(package_install_loc, .libPaths())
    .libPaths(paths)

}, error = function(e) report_error_and_exit("TMPDIR / library path setup", e))

# Load required libraries
tryCatch({
    message("Loading required libraries ...")
    suppressPackageStartupMessages({
        # Loaded first + patched before RiboseQC
        library(txdbmaker)
        library(GenomicFeatures)
    })

    ns <- asNamespace("GenomicFeatures")
    if (bindingIsLocked("makeTxDbFromGFF", ns)) {
        unlockBinding("makeTxDbFromGFF", ns)
    }
    assign("makeTxDbFromGFF", txdbmaker::makeTxDbFromGFF, envir = ns)
    lockBinding("makeTxDbFromGFF", ns)
    message("Patched GenomicFeatures::makeTxDbFromGFF -> txdbmaker::makeTxDbFromGFF")

    suppressPackageStartupMessages({
        library(RiboseQC)
        library(BSgenome)
        library(Biostrings)
    })

}, error = function(e) report_error_and_exit("library loading", e))

# Prepare annotation files
tryCatch({
    message("Preparing annotation ...")
    suppressWarnings(
        prepare_annotation_files(
            annotation_directory = savedir,
            twobit_file          = twobit_file,
            gtf_file             = gtf,
            genome_seq           = genome,
            annotation_name      = annot_name,
            forge_BSgenome       = TRUE
        )
    )
    message("Annotation preparation completed successfully")

}, error = function(e) report_error_and_exit("prepare_annotation_files()", e))

message("Finished creating index files")
quit(status = 0, save = "no")