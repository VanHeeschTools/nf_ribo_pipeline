#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(ORFquant)
  library(RiboseQC)
  library(rmarkdown)
})

# Obtain input arguments
args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
name <- args[2]
rannot <- args[3]
package_install_loc <- args[4]
readlength_choice_method <- args[5]


message("Running RiboseQC")
message(paste("Chosen readlength_choice_method: ", readlength_choice_method ))
paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Find Pandoc in container
rmarkdown::find_pandoc(dir = "/usr/src/pandoc/bin")

# Define functions
riboseqc_analysis <- function(bam, rannot, name) {
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
        readlength_choice_method = readlength_choice_method,
        report_file = name)
      message("Successfully executed RiboseQC")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  )
}
# Run script
riboseqc_analysis(bam, rannot, name)
