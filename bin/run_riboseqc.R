#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(ORFquant)
  library(RiboseQC)
  library(rmarkdown)
})

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
name <- args[2]
rannot <- args[3]
pandoc_dir <- args[4]
annotation_package <- args[5]
package_install_loc <- args[6]

paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Define variables --------------------------------------------------------
find_pandoc(dir = pandoc_dir)

# Define functions --------------------------------------------------------
riboseqc_analysis <- function(bam, rannot, name) {

  print(bam)

  tryCatch(

    expr = {

      RiboseQC_analysis(annotation_file = rannot,
                        bam_files = bam,
                        read_subset = FALSE,
                        dest_names = name,
                        rescue_all_rls = FALSE,
                        fast_mode = FALSE,
                        create_report = FALSE,
                        sample_names = NA,
                        report_file = name)

      message("Successfully executed the call.")

    },

    error = function(e){
      message('Caught an error!')
      print(e)
    }

  )
}

# See if script works without package
#if (!require(basename(annotation_package), character.only = TRUE)) {
#   install.packages(annotation_package,
#                    repos = NULL, type = "source")
# }

# Run script --------------------------------------------------------------
riboseqc_analysis(bam, rannot, name)
