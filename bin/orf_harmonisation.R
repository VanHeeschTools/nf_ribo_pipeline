#!/usr/bin/env Rscript

# =============================================================================
# 01 | LOAD LIBRARIES ----
# =============================================================================
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(rtracklayer)
})

# =============================================================================
# 02 | LOAD DATA ----
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
#orfquant_table <- args[1]
#price_table <- args [2]
# Gives NA value if RiboTIE input equals NULL
#ribotie_table <- ifelse(length(args) >= 3 && args[3] != "null", args[3], NA)
orfcaller_tables <- args

# =============================================================================
# 03 | FUNCTIONS ----
#   * Define all functions
# =============================================================================

#' Read annotated ORF tables
#'
#' This function takes a vector of file paths pointing to ORF quantification
#' tables (CSV files) reads each file as data.frame and puts it into the output list
#'
#' @param orfcaller_tables A character vector of file paths to ORF table files to be read.
#'
#' @return A list of loaded ORF table data frames
read_orf_tables <- function(orfcaller_tables) {
  
  orf_tables <- lapply(orfcaller_tables, function(f) {
    read.delim(
      f,
      sep = ",",
      colClasses = c(chrm = "character"),
      stringsAsFactors = FALSE
    )
  })
  
  return(orf_tables)
}

#' Filter identical ORFs preferring specific callers
#'
#' @param orfs Combined ORF data frame from all callers
#' @return Filtered ORF data frame with caller preference and count of callers per ORF
orf_filter <- function(orfs){
  caller_order <- c("ORFquant", "PRICE", "RiboTIE")
  
  filtered_table <- orfs %>%
    group_by(tx_id, protein_seq, starts, ends) %>%
    mutate(
      caller_count = n_distinct(orfcaller),  # Amount of callers ORF occurs in
      pref = match(orfcaller, caller_order), # Add preference column
    ) %>%
    slice_min(pref, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-pref) # Remove helper column
  
  return(filtered_table)
}

#' Sort ORF table by chromosome and coordinates
#'
#' @param filtered_orfs Data frame of filtered ORFs
#' @return Sorted data frame by chromosome and start/end positions
sort_orfs <- function(filtered_table){
  filtered_table_sorted <- filtered_table %>%
    # Arrange by chromosome order, then start and end
    dplyr::arrange(chrm, orf_start, orf_end) %>%

  return(filtered_table_sorted)
}

#' Obtain IDs of ORFs removed by filtering and write to file
#'
#' @param orfs Original combined ORF data frame
#' @param filtered_orfs Filtered ORF data frame
obtain_removed_orf_ids <- function(orfs, filtered_orfs){
  removed_orf_ids <- anti_join(orfs, filtered_orfs, by = "orf_id") %>%
    pull(orf_id)
  
  # Save removed orf_ids to a text file
  writeLines(removed_orf_ids, "removed_orf_ids.txt")
}

#' Write sorted data frame to CSV file
#'
#' @param sorted_df Sorted data frame of ORFs
#' @param output_file Output CSV file path
write_results <- function(sorted_df, output_file){
  write.table(sorted_df, file = output_file,
              sep = ",",
              quote = F,
              row.names = F)
}

# =============================================================================
# 04 | EXTRA OUTPUT FUNCTIONS ----
# =============================================================================

#' Write a protein FASTA file from a dataframe
#'
#' This function takes the orf_table dataframe and writes a compressed FASTA file (.fa.gz). 
#' Each `orf_id` is used as the FASTA header, and the corresponding `Protein` entry 
#' is written as the sequence, wrapped at 60 characters per line.
#'
#' @param sorted_df
#' @param fasta_file

write_orf_protein_fasta <- function(sorted_df, fasta_file) {
  # Open a gzipped file connection
  con <- gzfile(fasta_file, "w")
  on.exit(close(con))
  
  # Helper: wrap sequence into 60-char lines
  wrap_seq <- function(seq, width = 60) {
    paste(strwrap(seq, width = width), collapse = "\n")
  }
  
  # Build FASTA entries
  fasta_entries <- paste0(">", sorted_df$orf_id, "\n",
                          vapply(sorted_df$protein_seq, wrap_seq, character(1)))
  
  # Write to file
  writeLines(fasta_entries, con)
}

#' Write a DNA FASTA file from a dataframe
#'
#' This function takes the orf_table dataframe and writes a compressed FASTA file (.fa.gz). 
#' Each `orf_id` is used as the FASTA header, and the corresponding `DNA` entry 
#' is written as the sequence, wrapped at 60 characters per line.
#'
#' @param sorted_df
#' @param fasta_file

write_orf_dna_fasta <- function(sorted_df, fasta_file) {
  # Open a gzipped file connection
  con <- gzfile(fasta_file, "w")
  on.exit(close(con))
  
  # Helper: wrap sequence into 60-char lines
  wrap_seq <- function(seq, width = 60) {
    paste(strwrap(seq, width = width), collapse = "\n")
  }
  
  # Build FASTA entries
  fasta_entries <- paste0(">", sorted_df$orf_id, "\n",
                          vapply(sorted_df$dna_seq, wrap_seq, character(1)))
  
  # Write to file
  writeLines(fasta_entries, con)
}

#' Turn the sorted and filtered harmonised ORF table into a gtf-like file
#' 
#' @param sorted_df data.frame produced by orf_filter()
convert_to_gtf <- function(sorted_df, gtf_output_file) {

  # Handle transcript rows
  transcripts <- sorted_df %>%
    mutate(
      start = orf_start, 
      end = orf_end,
      feature    = "transcript",
      score      = ".",
      frame      = ".",
      attributes = paste0('transcript_id "', orf_id, '"; gene_id "', 
                          gene_id, '"; gene_name "', gene_name,
                          '"; gene_biotype "', gene_biotype, 
                          '"; ORF_category "', orf_biotype_single, 
                          '"; ORFcaller "', orfcaller, '";')
    ) %>%
    # Orf_id will be used to join with cds rows, and is removed afterwards
    dplyr::select(chrm, orfcaller, feature, start, end, 
                  score, strand, frame, attributes, orf_id) %>%
    arrange(chrm, start) # Sort based on genomic location
  
  cds <- sorted_df %>%
    mutate(
      starts = strsplit(as.character(starts), "_"),
      ends   = strsplit(as.character(ends), "_")
    ) %>%
    unnest(c(starts, ends)) %>%
    mutate(
      start = as.integer(starts),
      end   = as.integer(ends)
    ) %>%
    mutate(
      feature    = "CDS",
      score      = ".",
      frame      = ".",
      attributes = paste0('transcript_id "', orf_id, '"; gene_id "', 
                          gene_id, '"; gene_name "', gene_name, 
                          '"; gene_biotype "', gene_biotype,
                          '"; orf_biotype "', orf_biotype_single, 
                          '"; orfcaller "', orfcaller, '";')
    ) %>%
    dplyr::select(chrm, orfcaller, feature, start, end, 
                  score, strand, frame, attributes, orf_id)
  
  
  # Combine transcript rows with their corresponding CDS rows
  # The CDS rows are sorted by strand
  gtf_list <- lapply(1:nrow(transcripts), function(i) {
    tx <- transcripts[i, ] # current transcript row index
    tx_cds <- cds %>% filter(orf_id == tx$orf_id) # CDS rows for this transcript
    
    # Combine transcript and its CDS rows
    bind_rows(tx, tx_cds)
  })
  
  # Combine all transcript+CDS groups into a single tibble
  gtf_out <- do.call(rbind, gtf_list) %>%
    dplyr::select(-orf_id) %>% # remove helper column used for grouping
    mutate(chrm = as.character(chrm)) %>%  # convert factor back to character
    # Write output gtf file
    write.table(gtf_output_file, sep = "\t", quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
}

#' Generate MultiQC table of ORF categories per ORFcaller
#'
#' @param orf_dfs Named list of data frames, each corresponding to an ORFcaller
multiqc_orfcaller_table <- function(orf_dfs) {

  all_orf_types <- unique(unlist(lapply(orf_dfs, function(df) df$orf_biotype_single)))
  
  count_tables <- lapply(orf_dfs, function(df) {
    caller_name <- unique(df$orfcaller)
    
    counts <- df %>%
      count(orf_biotype_single) %>%
      complete(orf_biotype_single = all_orf_types, fill = list(n = 0)) %>%
      arrange(match(orf_biotype_single, all_orf_types)) %>%
      dplyr::select(n) %>%
      t() %>%
      as.data.frame()
    
    colnames(counts) <- all_orf_types
    dplyr::mutate(counts, ORFcaller = caller_name, .before = 1)
  })
  
  multiqc_table <- bind_rows(count_tables)
  sorted_cols <- names(sort(colMeans(multiqc_table[, -1]), decreasing = TRUE))
  multiqc_table %>%
    dplyr::select(ORFcaller, dplyr::all_of(sorted_cols))
  
  write.table(multiqc_table, "orfcaller_orf_categories_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Generate MultiQC table for merged ORF data
#'
#' @param sorted_df Sorted data frame of merged ORFs
multiqc_merged_table <- function(sorted_df){
  merged_orf_types <- unique(sorted_df$orf_biotype_single)
  
  merged_counts <- sorted_df %>%
    count(orf_biotype_single) %>%
    complete(orf_biotype_single = merged_orf_types, fill = list(n = 0)) %>%
    arrange(match(orf_biotype_single, merged_orf_types)) %>%
    dplyr::select(n) %>%
    t() %>%
    as.data.frame()
  
  colnames(merged_counts) <- merged_orf_types
  merged_counts <- mutate(merged_counts, ORFcaller = "Merged_ORFcallers", .before = 1)
  
  # Sort columns by descending average counts (which is just the counts here)
  sorted_cols_merged <- names(sort(colMeans(merged_counts[, -1]), decreasing = TRUE))
  merged_counts <- merged_counts %>% dplyr::select(ORFcaller, all_of(sorted_cols_merged))
  
  write.table(merged_counts, "merged_orf_categories_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Create MultiQC-style table of counts by caller_count (found in n callers)
#'
#' @param sorted_df data.frame produced by orf_filter(), must contain caller_count
#' @param outfile output filename (tab-separated). Default: "caller_count_mqc.txt"
#' @return data.frame showing in how many ORFcallers the ORF is found
multiqc_caller_count <- function(sorted_df, outfile = "caller_count_mqc.txt") {
  
  # Determine maximum number of callers represented
  max_calls <- max(sorted_df$caller_count, na.rm = TRUE)
  if (!is.finite(max_calls) || max_calls < 1) max_calls <- 1L
  
  # Count how many ORFs occur in exactly 1,2,.,n callers
  counts_vec <- vapply(1:max_calls, function(k) {
    sum(sorted_df$caller_count == k, na.rm = TRUE)
  }, integer(1))
  
  # Build a one-row data.frame with readable column names
  col_names <- paste0("found_in_", 1:max_calls)
  out_df <- as.data.frame(t(counts_vec), stringsAsFactors = FALSE)
  names(out_df) <- col_names
  out_df <- tibble::add_column(out_df, ORFcaller = "Merged_ORFcallers", .before = 1)
  
  # Write tab-separated file
  write.table(out_df, "merged_orf_caller_count_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

# =============================================================================
# 05 | RUN HARMONISATION ----
# =============================================================================

# Step 1: Load and merge ORFS

#orfquant_orfs <- read.delim(orfquant_table, sep = ",", colClasses = c(chrm = "character"), stringsAsFactors = FALSE)
#price_orfs <- read.delim(price_table, sep = ",",, colClasses = c(chrm = "character"), stringsAsFactors = FALSE)

#if (!is.na(ribotie_table)) {
#  ribotie_orfs <- read.delim(ribotie_table, sep = ",", colClasses = c(chrm = "character"), stringsAsFactors = FALSE)

  # Merge tables with RiboTIE results
#  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs, ribotie_orfs)
#} else {
  # Merge tables without RiboTIE results
#  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs)
#}

# Step 1: Load and merge ORFS
loaded_orf_tables <- read_orf_tables(orfcaller_tables)
orfs <- dplyr::bind_rows(loaded_orf_tables)

# Step 2: Remove identical ORFs between ORFcallers
filtered_table <- orf_filter(orfs)

# Step 3: Sorted ORF table
sorted_df <- sort_orfs(filtered_table)

# Step 4: Obtain removed duplicates ORF ids and write IDs to txt file
removed_orfs <- obtain_removed_orf_ids(orfs, filtered_table)

# Step 5: Write ORF protein and DNA sequences to fasta files
write_orf_protein_fasta(sorted_df, "orf_protein_sequences.fa.gz")
write_orf_dna_fasta(sorted_df, "orf_dna_sequences.fa.gz")

# Remove DNA-seq from harmonised ORF table
sorted_df <- sorted_df %>%
  dplyr::select(-dna_seq)

# Step 6: Write harmonised ORF table to csv file
write_results(sorted_df, "harmonised_orf_table.csv")

# Step 7: Convert harmonised ORF table to gtf file
convert_to_gtf(sorted_df, "harmonised_orf_table.gtf")

# For testing purposes sort and write unfiltered harmonised orf table
sorted_unfiltered_df <- sort_orfs(orfs) %>%
  dplyr::select(-dna_seq)

write_results(sorted_unfiltered_df, "unfiltered_harmonised_orf_table.csv")

# =============================================================================
# 06 | CREATE MULTIQC TABLES ----
# =============================================================================

if (!is.na(ribotie_table)) {
  orf_dfs <- list(
    ORFquant = orfquant_orfs,
    PRICE    = price_orfs,
    RiboTIE  = ribotie_orfs
  )
} else {
  orf_dfs <- list(
    ORFquant = orfquant_orfs,
    PRICE    = price_orfs
  )
}

multiqc_orfcaller_table(loaded_orf_tables)
multiqc_merged_table(sorted_df)
multiqc_caller_count(sorted_df)