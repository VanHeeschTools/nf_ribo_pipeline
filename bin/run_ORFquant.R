#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(ORFquant)
  library(RiboseQC)
  library(rmarkdown)
})

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
for_ORFquant_file <- args[1]
rannot <- args[2]
cpu <- args[3]
pandoc_dir <- args[4]
annotation_package <- args[5]
package_install_loc <- args[6]


# REMOVE AFTER TESTING
source("/hpc/pmc_vanheesch/projects/evanderwerf/nextflow_pipelines/riboseq_pipeline/orfquant_script/orfquant.R")
source("/hpc/pmc_vanheesch/projects/evanderwerf/nextflow_pipelines/riboseq_pipeline/orfquant_script/FaFile_Circ.R")

pool_id = "output"
paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Define variables --------------------------------------------------------
ORFquant_output_file <- paste0(pool_id, "_final_ORFquant_results")
ORFquant_plot_data <- paste0(ORFquant_output_file, "_plots/",
                             pool_id, "_ORFquant_plots_RData")
find_pandoc(dir = pandoc_dir)

# Define functions --------------------------------------------------------
ORFquant_analysis <- function(for_ORFquant_file,
                              sample_id,
                              ORFquant_output_file,
                              ORFquant_plot_data,
                              rannot,
                              cpu,
                              html_report) {

  gene_names <- NA
  gene_ids <- NA
  
  if (!file.exists(ORFquant_output_file)) {

    run_ORFquant(for_ORFquant_file = for_ORFquant_file,
      annotation_file = rannot,
      n_cores = cpu,
      prefix = pool_id,
      gene_name = gene_names,
      gene_id = gene_ids,
      genomic_region = NA,
      write_temp_files = TRUE,
      write_GTF_file = TRUE,
      write_protein_fasta = TRUE,
      interactive = FALSE,
      stn.orf_find.all_starts = TRUE,
      stn.orf_find.nostarts = FALSE,
      stn.orf_find.start_sel_cutoff = NA,
      stn.orf_find.start_sel_cutoff_ave = 0.5,
      stn.orf_find.cutoff_fr_ave = 0.5,
      stn.orf_quant.cutoff_cums = NA,
      unique_reads_only = FALSE,
      canonical_start_only = TRUE,
      stn.orf_quant.scaling = "total_Psites"
    )
  } else {
    print("ORFquant file already exists")
  }

  plot_ORFquant_results(for_ORFquant_file = for_ORFquant_file,
    ORFquant_output_file = ORFquant_output_file,
    annotation_file = rannot
  )

}

# Run script --------------------------------------------------------------
# See if script works without annotation
# if (!require(basename(annotation_package), character.only = TRUE)) {
#   install.packages(annotation_package,
#                    repos = NULL, type = "source")
# }

ORFquant_analysis(for_ORFquant_file = for_ORFquant_file,
                  sample_id = pool_id,
                  ORFquant_output_file = ORFquant_output_file,
                  ORFquant_plot_data = ORFquant_plot_data,
                  rannot = rannot,
                  cpu = cpu,
                  html_report = html_report)