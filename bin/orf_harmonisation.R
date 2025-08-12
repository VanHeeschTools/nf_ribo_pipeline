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
orfquant_table <- args[1]
price_table <- args [2]
# Gives NA value if RiboTIE input equals NULL
ribotie_table <- ifelse(length(args) >= 3 && args[3] != "null", args[3], NA)

# =============================================================================
# 03 | FUNCTIONS ----
#   * Define all functions
# =============================================================================

#' Read ORFquant CSV and preprocess
#'
#' @param orfquant_table File path to ORFquant CSV file
#' @return Data frame with ORFquant ORFs
read_orfquant_csv <- function(orfquant_table){
  orfquant_orfs <- read.delim(orfquant_table, sep = ",") %>%
    # Remove unnecessary columns
    dplyr::select(-c("transcript_biotype")) %>%
    # Add orf_caller column for identification
    dplyr::mutate(orf_caller = "ORFquant")

  return(orfquant_orfs)
}

#' Read Price CSV and preprocess
#'
#' @param price_table File path to Price CSV file
#' @return Data frame with Price ORFs
read_price_csv <- function(price_table){
  price_orfs <- read.delim(price_table, sep = ",") %>%
    # Remove unnecessary columns
    dplyr::select(-c("seqnames", "start", "end", "width")) %>%
    # Add orf_caller column for identification
    dplyr::mutate(orf_caller = "PRICE") 
  
  return(price_orfs)
}

#' Read RiboTIE CSV and preprocess
#'
#' @param ribotie_table File path to RiboTIE CSV file
#' @return Data frame with RiboTIE ORFs
read_ribotie_csv <- function(ribotie_table){
  ribotie_orfs <- read.delim(ribotie_table, sep = ",") %>%
    # Remove unnecessary columns
    dplyr::select(-c("seqname",
                    "CDS_coords",
                    "ORF_type",
                    "ORF_len")) %>%
    # Add orf_caller column for identification
    dplyr::mutate(orf_caller = "RiboTIE")
  
  return(ribotie_orfs)
}

#' Filter identical ORFs preferring specific callers
#'
#' @param orfs Combined ORF data frame from all callers
#' @return Filtered ORF data frame with caller preference and count of callers per ORF
orf_filter <- function(orfs){
  caller_order <- c("ORFquant", "PRICE", "RiboTIE")
  
  filtered_table <- orfs %>%
    group_by(gene_id, Protein, ORF_ranges) %>%
    mutate(
      caller_count = n_distinct(orf_caller), # Amount of callers ORF occurs in
      pref = match(orf_caller, caller_order), # Add preference column
    ) %>%
    slice_min(pref, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-pref) # Remove helper column
  
  return(filtered_table)
}

#filtered_table <- orfs %>%
#  group_by(gene_id, Protein, ORF_ranges) %>%
#  slice_sample(n = 1) %>% 
#  ungroup()

#' Sort ORF table by chromosome and coordinates
#'
#' @param filtered_orfs Data frame of filtered ORFs
#' @return Sorted data frame by chromosome and start/end positions
sort_orfs <- function(filtered_table){
filtered_table_sorted <- filtered_table %>%
  # Split ORF_ranges into chr and coordinate range columns
  tidyr::separate(ORF_ranges, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  dplyr::mutate(
      # Clean chromosome string and convert to numeric when possible
      chr_clean = gsub("^chr", "", chr),
      chr_numeric = suppressWarnings(as.numeric(chr_clean)),
      chr_order = ifelse(is.na(chr_numeric), chr_clean, chr_numeric)
  ) %>%
  # Arrange by chromosome order, then start and end
  dplyr::arrange(chr_order, start, end) %>%
  # Remove helper columns
  dplyr::select(-chr_clean, -chr_numeric, -chr_order, -start, -end, -chr)

  return(filtered_table_sorted)
}

#' Obtain IDs of ORFs removed by filtering and write to file
#'
#' @param orfs Original combined ORF data frame
#' @param filtered_orfs Filtered ORF data frame
obtain_removed_orf_ids <- function(orfs, filtered_orfs){
  removed_orf_ids <- anti_join(orfs, filtered_table, by = c("gene_id", "gene_name", "Protein", "transcript_id", "orf_id")) %>%
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

#' Generate MultiQC table of ORF categories per ORFcaller
#'
#' @param orf_dfs Named list of data frames, each corresponding to an ORFcaller
multiqc_orfcaller_table <- function(orf_dfs) {
  all_orf_types <- unique(unlist(lapply(orf_dfs, \(df) unique(df$orf_category_new))))

  count_tables <- lapply(names(orf_dfs), function(caller_name) {
    df <- orf_dfs[[caller_name]]
    counts <- df %>%
      count(orf_category_new) %>%
      complete(orf_category_new = all_orf_types, fill = list(n = 0)) %>%
      arrange(match(orf_category_new, all_orf_types)) %>%
      select(n) %>%
      t() %>%
      as.data.frame()
    colnames(counts) <- all_orf_types
    mutate(counts, ORFcaller = caller_name, .before = 1)
  })

  multiqc_table <- bind_rows(count_tables)
  sorted_cols <- names(sort(colMeans(multiqc_table[ , -1]), decreasing = TRUE))
  multiqc_table <- multiqc_table %>% select(ORFcaller, all_of(sorted_cols))
  
  write.table(multiqc_table, "orfcaller_orf_categories_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Generate MultiQC table for merged ORF data
#'
#' @param sorted_df Sorted data frame of merged ORFs
multiqc_merged_table <- function(sorted_df){
  merged_orf_types <- unique(sorted_df$orf_category_new)

  merged_counts <- sorted_df %>%
    count(orf_category_new) %>%
    complete(orf_category_new = merged_orf_types, fill = list(n = 0)) %>%
    arrange(match(orf_category_new, merged_orf_types)) %>%
    select(n) %>%
    t() %>%
    as.data.frame()

  colnames(merged_counts) <- merged_orf_types
  merged_counts <- mutate(merged_counts, ORFcaller = "Merged_ORFcallers", .before = 1)

  # Sort columns by descending average counts (which is just the counts here)
  sorted_cols_merged <- names(sort(colMeans(merged_counts[, -1]), decreasing = TRUE))
  merged_counts <- merged_counts %>% select(ORFcaller, all_of(sorted_cols_merged))

  write.table(merged_counts, "merged_orf_categories_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Create MultiQC-style table of counts by caller_count (found in 1,2,... callers)
#'
#' @param sorted_df data.frame produced by orf_filter(), must contain caller_count
#' @param outfile output filename (tab-separated). Default: "caller_count_mqc.txt"
#' @return data.frame (one row) with columns found_in_1, found_in_2, ...
multiqc_caller_count <- function(sorted_df, outfile = "caller_count_mqc.txt") {

  # determine maximum number of callers represented
  max_calls <- max(sorted_df$caller_count, na.rm = TRUE)
  if (!is.finite(max_calls) || max_calls < 1) max_calls <- 1L

  # count how many ORFs occur in exactly 1,2,...,max_calls callers
  counts_vec <- vapply(1:max_calls, function(k) {
    sum(sorted_df$caller_count == k, na.rm = TRUE)
  }, integer(1))

  # build a one-row data.frame with readable column names
  col_names <- paste0("found_in_", 1:max_calls)
  out_df <- as.data.frame(t(counts_vec), stringsAsFactors = FALSE)
  names(out_df) <- col_names
  out_df <- tibble::add_column(out_df, ORFcaller = "Merged_ORFcallers", .before = 1)

  # write tab-separated file
  write.table(out_df, "merged_orf_caller_count_mqc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
}

# =============================================================================
# 04 | RUN HARMONISATION ----
# =============================================================================

# Step 1: Load and merge ORFS
orfquant_orfs <- read_orfquant_csv(orfquant_table)
price_orfs <- read_price_csv(price_table)
if (!is.na(ribotie_table)) {
  ribotie_orfs <- read_ribotie_csv(ribotie_table)
  # Merge tables with RiboTIE results
  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs, ribotie_orfs)
} else {
  # Merge tables without RiboTIE results
  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs)
}

# Step 2: Remove identical ORFs between ORFcallers
filtered_table <- orf_filter(orfs)

# Step 3: Sorted ORF table
sorted_df <- sort_orfs(filtered_table)

# Step 4: Obtain removed duplicates ORF ids and write IDs to txt file
removed_orfs <- obtain_removed_orf_ids(orfs, filtered_table)

# Step 5: Write harmonised ORF table to csv file
write_results(sorted_df, "harmonised_table.csv")

# For testing purposes sort and write unfiltered harmonised orf table
sorted_unfiltered_df <- sort_orfs(orfs)
write_results(sorted_unfiltered_df, "unfiltered_harmonised_table.csv")

# =============================================================================
# 05 | CREATE MULTIQC TABLES ----
# =============================================================================

if (!is.na(ribotie_table)) {
  orf_dfs <- list(
    ORFquant = orfquant_orfs,
    PRICE    = price_orfs,
    RiboTIE  = ribotie_orfs
  )
}else{
  orf_dfs <- list(
    ORFquant = orfquant_orfs,
    PRICE    = price_orfs
  )
}
multiqc_orfcaller_table(orf_dfs)
multiqc_merged_table(sorted_df)
multiqc_caller_count(sorted_df)
