#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(stringr)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ------------------------------------------------------------------------------
# Load RiboseQC output files and initialize color palette
# ------------------------------------------------------------------------------
input_files <- commandArgs(trailingOnly = TRUE)

if (length(input_files) == 0) stop("No input files provided.")

#riboseqc_files <- list.files("riboseqc", pattern = "results_RiboseQC_all", full.names = T, recursive = T)

# ------------------------------------------------------------------------------
# Initialize empty data frames for collecting statistics from each sample
# ------------------------------------------------------------------------------
summary_P_sites_df <- data.frame()
summary_reads_df <- data.frame()
inframe_df <- data.frame()
read_cats_df <- data.frame()
cds_reads_df <- data.frame()

# ------------------------------------------------------------------------------
# Process each RiboseQC results file
# ------------------------------------------------------------------------------
for (fname in input_files) {
  parts <- str_split(basename(fname), "_")[[1]]
  sample_id <- paste(parts[1:(length(parts) - 3)], collapse = "_")
   
  message("Loading ", sample_id)
  load(fname)
  
  # Extract data from loaded RiboseQC result
  summary_P_sites_sample <- as.data.frame(res_all$summary_P_sites)
  summary_P_sites_sample$sample_id <- sample_id
  
  summary_reads_sample <- as.data.frame(t(colSums(as.data.frame(res_all$read_stats$reads_summary_unq$nucl))))
  rownames(summary_reads_sample) <- sample_id
  
  inframe_sample <- as.data.frame(t(res_all$selection_cutoffs$analysis_frame_cutoff$nucl$all$frames_res))
  rownames(inframe_sample) <- sample_id
  
  read_cats_sample <- as.data.frame(t(rowSums(as.data.frame(res_all$read_stats$reads_summary$nucl))))
  rownames(read_cats_sample) <- sample_id
  
  cds_reads_sample <- data.frame(reads = sum(res_all$read_stats$counts_cds_genes_unq$reads))
  rownames(cds_reads_sample) <- sample_id
  
  # Merge into master data frames
  summary_P_sites_df <- rbind(summary_P_sites_df, summary_P_sites_sample)
  summary_reads_df <- bind_rows(summary_reads_df, summary_reads_sample)
  inframe_df <- rbind(inframe_df, inframe_sample)
  read_cats_df <- rbind(read_cats_df, read_cats_sample)
  cds_reads_df <- rbind(cds_reads_df, cds_reads_sample)
}

# ------------------------------------------------------------------------------
# Frame Preference Table (29-nt reads, nucleotide composition)
# Output: MultiQC-compatible table
# ------------------------------------------------------------------------------
summary_P_sites_df_s <- summary_P_sites_df  # Ensure scoped correctly
frame_29 <- summary_P_sites_df_s %>%
  filter(read_length == 29, comp == "nucl") %>%
  select(sample_id, frame_preference)

writeLines(c(
  "Sample\tFramePreference",
  paste(frame_29$sample_id, round(frame_29$frame_preference, 2), sep = "\t")
), con = "riboseqc_frame_29nt_mqc.txt")

# ------------------------------------------------------------------------------
# Read Category Counts Table
# Output: MultiQC-compatible table
# ------------------------------------------------------------------------------
# Prepare final table for MultiQC output
read_cats_raw <- read_cats_df
read_cats_raw$Sample <- rownames(read_cats_raw)
read_cats_raw <- read_cats_raw[, c(ncol(read_cats_raw), 1:(ncol(read_cats_raw) - 1))]

header <- paste(colnames(read_cats_raw), collapse = "\t")
lines <- apply(read_cats_raw, 1, function(row) paste(row, collapse = "\t"))

writeLines(c(
  header,
  lines
), con = "riboseqc_read_categories_counts_mqc.txt")
