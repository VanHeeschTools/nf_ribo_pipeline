#!/usr/bin/env Rscript

#' Ribo-seq GTF processing and writing
#'
#' Load ORF tables and caller GTFs, filter by ORF IDs, merge, and write GTFs
#' with optional transcript rows.

suppressPackageStartupMessages({
    library(rtracklayer)
    library(dplyr)
    library(tibble)
    library(purrr)
    library(stringr)
})

# ----------------------
# Load input files
# ----------------------

args <- commandArgs(trailingOnly = TRUE)
orf_table_csv <- args[1]
orfquant_gtf <- args[2]
price_gtf <- args [3]
# Gives NA value if RiboTIE input equals NULL
ribotie_gtf <- ifelse(length(args) >= 4, args[4], NA)

# Load ORF_table
orf_table <- read.csv(orf_table_csv)

# ----------------------
# Functions
# ----------------------

#' Filter GTF rows by ORF IDs from a given caller
#' 
#' @param gtf_file Path to GTF file
#' @param orf_list_full Data frame of ORF IDs and callers
#' @param orfcaller_name Name of the ORF caller
#' @return tibble of filtered GTF rows
obtain_orfcaller_gtf <- function(gtf_file, orf_list_full, orfcaller_name) {
  orf_ids <- orf_list_full %>% 
    dplyr::filter(orfcaller == orfcaller_name) %>% 
    dplyr::select(orf_id)
  
  rtracklayer::import(gtf_file) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(-dplyr::any_of("start_codon")) %>% 
    semi_join(orf_ids, by = c("ORF_id" = "orf_id"))
}

#' Alter GTF dataframe
#' 
#' @param gtf_df GTF tibble
#' @param orf_table ORF table data frame
#' @return Altered GTF tibble
alter_gtf <- function(gtf_df, orf_table) {

  gtf_df <- gtf_df %>% 
    left_join(orf_table %>% 
                dplyr::select(orf_id, orf_biotype_single),
                by = c("ORF_id" = "orf_id"))

  return(gtf_df)
}

#' Write GTF with transcript rows which contain ORF_id as transcript and the CDS rows
#' 
#' @param df GTF tibble
#' @param file Output file path
write_gtf_with_transcripts <- function(df, file) {
  
  df <- merged_gtf
  
  stopifnot(all(c("seqnames","start","end","strand","ORF_id",
                  "gene_id","gene_name","orf_biotype_single","source") %in% names(df)))
  
  # Determine chromosome order
  seq_levels <- sort(unique(df$seqnames))
  
  # Handle transcript rows
  transcripts <- df %>%
    group_by(ORF_id, seqnames, strand, gene_id, gene_name, orf_biotype_single, gene_biotype, source) %>%
    summarise(start = min(start), end = max(end), .groups = "drop") %>%
    mutate(
      feature    = "transcript",
      score      = ".",
      frame      = ".",
      attributes = paste0('transcript_id "', ORF_id, '"; gene_id "', gene_id, '"; gene_name "', gene_name,'"; gene_biotype "', gene_biotype,
                          '"; ORF_category "', orf_biotype_single, '"; ORFcaller "', source, '";'),
      seqnames   = factor(seqnames, levels = seq_levels)
    ) %>%
    dplyr::select(seqnames, source, feature, start, end, score, strand, frame, attributes, ORF_id) %>%
    arrange(seqnames, start) # Sort based on genomic location
  
  # Handle CDS rows
  cds <- df %>%
    mutate(
      feature    = "CDS",
      score      = ".",
      frame      = ".",
      attributes = paste0('transcript_id "', ORF_id, '"; gene_id "', gene_id, '"; gene_name "', gene_name, '"; gene_biotype "', gene_biotype,
                          '"; ORF_category "', orf_biotype_single, '"; ORFcaller "', source, '";')
    ) %>%
    dplyr::select(seqnames, source, feature, start, end, score, strand, frame, attributes, ORF_id)
  
  # Combine transcript rows with their corresponding CDS rows
  # The CDS rows are sorted by strand
  gtf_list <- lapply(1:nrow(transcripts), function(i) {
    tx <- transcripts[i, ] # current transcript row index
    tx_cds <- cds %>% filter(ORF_id == tx$ORF_id)  # CDS rows for this transcript
    
    # Sort CDS according to strand direction
    tx_cds <- if (tx$strand == "+") arrange(tx_cds, start) else arrange(tx_cds, desc(start))
    
    # Combine transcript and its CDS rows
    bind_rows(tx, tx_cds)
  })
  
  # Combine all transcript+CDS groups into a single tibble
  gtf_out <- do.call(rbind, gtf_list) %>%
    dplyr::select(-ORF_id) %>% # remove helper column used for grouping
    mutate(seqnames = as.character(seqnames)) %>%  # convert factor back to character
    # Write output gtf file
    write.table(file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}

# ----------------------
# Run workflow
# ----------------------

# Load ORF_table csv file
orf_list_full <- orf_table %>% dplyr::select(orf_id, orfcaller)

# Obtain gtf rows of ORFquant ORFs and add/remove attributes
gtf_orfquant <- obtain_orfcaller_gtf(orfquant_gtf, orf_list_full, "ORFquant")
gtf_orfquant_altered <- alter_gtf(gtf_orfquant, orf_table)

# Obtain gtf rows of PRICE ORFs and add/remove attributes
gtf_price <- obtain_orfcaller_gtf(price_gtf, orf_list_full, "PRICE")
gtf_price_altered    <- alter_gtf(gtf_price, orf_table)

# Only handle RiboTIE is input files is given
if (!is.na(ribotie_gtf)) {
  gtf_ribotie  <- obtain_orfcaller_gtf(ribotie_gtf,  orf_list_full, "RiboTIE")
  gtf_ribotie_altered  <- alter_gtf(gtf_ribotie,  orf_table)
  
  # Merge all gtfs including RiboTIE
  merged_gtf <- bind_rows(gtf_orfquant_altered, gtf_price_altered, gtf_ribotie_altered) %>%
    arrange(seqnames, start, end, strand)
  
}else{
  # Merge Orfquant and PRICE gtfs
  merged_gtf <- bind_rows(gtf_orfquant_altered, gtf_price_altered) %>%
    arrange(seqnames, start, end, strand)
}

# Write gtf files
write_gtf_with_transcripts(merged_gtf, "combined_orfcallers.gtf") 

