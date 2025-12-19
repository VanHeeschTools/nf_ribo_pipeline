#!/usr/bin/env Rscript

message("Loading required libraries")
suppressPackageStartupMessages({
    library(RiboseQC)
    library(rmarkdown)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[-length(args)]
rmd_path <- args[length(args)]

input_sample_names <- gsub("_results_RiboseQC_all", "", basename(input_files))

# Find Pandoc in container
rmarkdown::find_pandoc(dir = "/usr/src/pandoc/bin")

# Paths
workdir = getwd()
tmp_dir   <- paste0(file.path(workdir, "tmp"), "/")
plots_dir <- paste0(file.path(workdir, "plots"), "/")
rds_dir   <- paste0(plots_dir, "rds/")
pdf_dir   <- paste0(plots_dir, "pdf/")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, showWarnings = FALSE)

sink("RiboseQC_report_text_output.txt")

message("Rendering HTML report")

knitr::knit_meta(class = NULL, clean = TRUE)

suppressWarnings(
    render(
        rmd_path,
        params = list(
            input_files = input_files,
            input_sample_names = input_sample_names,
            output_fig_path = plots_dir
        ),
        output_file = file.path(workdir, "RiboseQC_report.html"),
        intermediates_dir = file.path(workdir, "tmp"),
        knit_root_dir = getwd()
    )
)

gc()
sink()
