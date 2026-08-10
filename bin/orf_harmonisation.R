#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(rtracklayer)
})

# Load input data
args <- commandArgs(trailingOnly = TRUE)
orfcaller_tables <- args

# Define functions

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
      colClasses = c(chr = "character", orf_start = "integer", orf_end = "integer"),
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
  
  # Annotate if ORFs are found in each ORFcaller
  annotated <- orfs %>%
    group_by(tx_id, protein_seq, starts, ends) %>%
    mutate(
      pref = match(orfcaller, caller_order),
      found_in_ORFquant = any(orfcaller == "ORFquant"),
      found_in_PRICE    = any(orfcaller == "PRICE"),
      found_in_RiboTIE  = any(orfcaller == "RiboTIE"),
      caller_count = n_distinct(orfcaller)
    ) %>%
    ungroup()
  
  # Keep preferred ORF per group
  filtered_table <- annotated %>%
    group_by(tx_id, protein_seq, starts, ends) %>%
    slice_min(pref, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-pref)  # remove helper column
  
  # Drop found_in columns if all values are FALSE
  found_cols <- c("found_in_ORFquant", "found_in_PRICE", "found_in_RiboTIE")
  unused_cols <- found_cols[sapply(filtered_table[found_cols], function(x) all(!x))]
  filtered_table <- filtered_table %>% dplyr::select(-dplyr::all_of(unused_cols))
  
  return(filtered_table)
}

#' Sort ORF table by chromosome and coordinates
#'
#' @param filtered_orfs Data frame of filtered ORFs
#' @return Sorted data frame by chromosome and start/end positions
sort_orfs <- function(filtered_table){
  filtered_table_sorted <- filtered_table %>%
    # Arrange by chromosome order, then start and end
    dplyr::arrange(chr, orf_start, orf_end) %>%
    dplyr::select(-orfcaller)

  return(filtered_table_sorted)
}

#' Obtain IDs of ORFs removed by filtering and write to file
#'
#' @param orfs Original combined ORF data frame
#' @param filtered_orfs Filtered ORF data frame
obtain_removed_orf_ids <- function(orfs, filtered_orfs){
  removed_orf_ids <- anti_join(orfs, filtered_orfs, by = "orfcaller_orf_id") %>%
    pull(orfcaller_orf_id)
  
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

#' Write a protein FASTA file from a dataframe
#'
#' This function takes the orf_table dataframe and writes a compressed FASTA file (.fa.gz). 
#' Each `summary_id` is used as the FASTA header, and the corresponding `Protein` entry 
#' is written as the sequence, wrapped at 60 characters per line.
#'
#' @param sorted_df
#' @param fasta_file

write_orf_protein_fasta <- function(sorted_df, fasta_file, mstart = FALSE) {
    # Open a gzipped file connection
    con <- gzfile(fasta_file, "w")
    on.exit(close(con))
    
    # Helper: wrap sequence into 60-char lines
    wrap_seq <- function(seq, width = 60) {
      paste(strwrap(seq, width = width), collapse = "\n")
    }
    
    # Replace the first amino acid with M if mstart is true
    if (mstart) {
      # Substitute the first character of each protein sequence with "M"
      seq_df <- sorted_df %>%
        dplyr::mutate(protein_seq_out = sub("^.", "M", protein_seq))
    } else {
      # Keep protein sequence unchanged
      seq_df <- sorted_df %>%
        dplyr::mutate(protein_seq_out = protein_seq)
    }
    
    # Wrap each protein sequence into 60-char lines, one row at a time
    seq_df <- seq_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate(protein_seq_wrapped = wrap_seq(protein_seq_out)) %>%
      dplyr::ungroup()
    
    # Build FASTA entries: ">summary_id" header line followed by wrapped sequence
    fasta_entries <- paste0(">", seq_df$summary_id, "\n", seq_df$protein_seq_wrapped)
    
    # Write all entries to the gzipped fasta file
    writeLines(fasta_entries, con)
}


#' Write a DNA FASTA file from a dataframe
#'
#' This function takes the orf_table dataframe and writes a compressed FASTA file (.fa.gz). 
#' Each `summary_id` is used as the FASTA header, and the corresponding `DNA` entry 
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
  fasta_entries <- paste0(">", sorted_df$summary_id, "\n",
                          vapply(sorted_df$dna_seq, wrap_seq, character(1)))
  
  # Write to file
  writeLines(fasta_entries, con)
}

#' Turn the harmonised ORF table into a gtf-like file
#' 
#' @param sorted_df data.frame produced by orf_filter()
convert_to_gtf <- function(sorted_df, gtf_output_file) {

    # Find all columns that show if ORF is found in caller
    found_cols <- grep("^found_in_", names(sorted_df), value = TRUE)
    found_attrs <- Reduce(paste0, lapply(found_cols, function(col) {
      paste0(col, ' "', sorted_df[[col]], '"; ')
    }))

    # Handle transcript rows
    transcripts <- sorted_df %>%
    mutate(
        source = "Harmonised",
        start = as.integer(orf_start), 
        end   = as.integer(orf_end),
        feature    = "transcript",
        score      = ".",
        frame      = ".",
        attributes = paste0('transcript_id "', orf_id, '"; gene_id "', 
                            gene_id, '"; gene_name "', gene_name, 
                            '"; gene_biotype "', gene_biotype,
                            '"; ORF_id "', orf_id,
                            '"; ORF_biotype_single "', orf_biotype_single, '"; ',
                            found_attrs)
    ) %>%
    # Orf_id will be used to join with cds rows, and is removed afterwards
    dplyr::select(chr, source, feature, start, end, 
                    score, strand, frame, attributes, orf_id) %>%
    arrange(chr, start) # Sort based on genomic location

    cds <- sorted_df %>%
    mutate(
        source = "Harmonised",
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
                            '"; ORF_id "', orf_id,
                            '"; ORF_biotype_single "', orf_biotype_single, '"; ', 
                            found_attrs)
    ) %>%
    dplyr::select(chr, source, feature, start, end, 
                    score, strand, frame, attributes, orf_id)

    # Combine transcript rows with their corresponding CDS rows
    # The CDS rows are sorted by strand
    n_tx <- nrow(transcripts) 
    cds_by_orf <- split(cds, cds$orf_id)

    gtf_list <- lapply(1:n_tx, function(i) {
        tx <- transcripts[i, ]
        tx_cds <- cds_by_orf[[tx$orf_id]]
        if (is.null(tx_cds)) tx else bind_rows(tx, tx_cds)
    })

    # Combine all transcript+CDS groups into a single tibble
    gtf_out <- do.call(rbind, gtf_list) %>%
    dplyr::select(-orf_id) %>% # remove helper column used for grouping
    mutate(chr = as.character(chr)) # convert factor back to character

    # Write output gtf file
    write.table(gtf_out, file = gtf_output_file, sep = "\t", quote = FALSE, 
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

# Run functions

# Load  ORF tables
loaded_orf_tables <- read_orf_tables(orfcaller_tables)

# If more than one ORF table is found merge them and remove identical ORFs
if (length(orfcaller_tables) > 1) {
  # Merge ORF tables
  orfs <- dplyr::bind_rows(loaded_orf_tables)
  # Remove identical ORFs between ORFcallers
  filtered_table <- orf_filter(orfs)
  # Obtain removed duplicates ORF ids and write IDs to txt file
  obtain_removed_orf_ids(orfs, filtered_table)
} else {
  filtered_table <- loaded_orf_tables[[1]]
  # Obtain removed duplicates ORF ids and write IDs to txt file
  obtain_removed_orf_ids(filtered_table, filtered_table)
}

sorted_df <- sort_orfs(filtered_table)

# Write ORF protein and DNA sequences to fasta files
write_orf_protein_fasta(sorted_df, "orf_protein_sequences.fa.gz", FALSE)
write_orf_protein_fasta(sorted_df, "orf_protein_sequences_M_start.fa.gz", TRUE)
write_orf_dna_fasta(sorted_df, "orf_dna_sequences.fa.gz")

# Remove extra columns from harmonised ORF table
sorted_df <- sorted_df %>%
  dplyr::select(-dna_seq)

# Write harmonised ORF table to csv file
write_results(sorted_df, "harmonised_orf_table.csv")

# Convert harmonised ORF table to gtf file
convert_to_gtf(sorted_df, "harmonised_orf_table.gtf")

# Create MultiQC tables
multiqc_orfcaller_table(loaded_orf_tables)
multiqc_merged_table(sorted_df)
if (length(orfcaller_tables) > 1) {
  multiqc_caller_count(sorted_df)
}