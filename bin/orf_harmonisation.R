#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)
orfquant_table <- args[1]
price_table <- args [2]
# Gives NA value if RiboTIE input equals NULL
ribotie_table <- ifelse(length(args) >= 3 && args[3] != "null", args[3], NA)

# Load ORF annotation table from PRICE and ORFquant
# TODO: Need to double check for input consistency between old and new code
price_orfs <- read.delim(price_table, sep = ",") %>%
  # Create a new column 'ORF_ranges' that combines genomic coordinates
  # Remove unnecessary columns
  dplyr::select(-c("seqnames", "start", "end", "width")) %>%
  # Add a orf_caller column for better tracking
  dplyr::mutate(orf_caller = "PRICE") 

orfquant_orfs <- read.delim(orfquant_table, sep = ",") %>%
  # Remove unnecessary columns
  dplyr::select(-c("transcript_biotype")) %>%
  dplyr::mutate(orf_caller = "ORFquant")

# Load RiboTIE results if RiboTIE was run
if (!is.na(ribotie_table)) {
  ribotie_orfs <- read.delim(ribotie_table, sep = ",") %>%
    # Remove unnecessary columns
    dplyr::select(-c("transcript_len",
                    "seqname",
                    "CDS_has_annotated_start_codon",
                    "CDS_has_annotated_stop_codon",
                    "CDS_coords",
                    "ORF_type",
                    "stop_codon",
                    "ORF_len")) %>%
    dplyr::mutate(orf_caller = "RiboTIE")

  # Combine all ORF tables into a single dataset
  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs, ribotie_orfs)

} else {
  # Combine both ORF tables into a single dataset
  orfs <- dplyr::bind_rows(orfquant_orfs, price_orfs)

}

# PRICE and ORFquant can predict identical ORF sequences in the same range
# Identical ORFs should be filtered based on:
# Gene id, Identical protein sequence, ORF_ranges

filtered_table <- orfs %>%
  # Preference = 0 for PRICE, 1 for everything else
  mutate(pref = ifelse(orf_caller == "PRICE", 0L, 1L)) %>%
  group_by(gene_id, Protein, ORF_ranges) %>%
  slice_min(pref, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-pref)

# Sort the table based on chomosome and coords
filtered_table_sorted <- filtered_table %>%
  tidyr::separate(ORF_ranges, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  dplyr::mutate(
      chr_clean = gsub("^chr", "", chr),
      chr_numeric = suppressWarnings(as.numeric(chr_clean)),
      chr_order = ifelse(is.na(chr_numeric), chr_clean, chr_numeric)
  ) %>%
  dplyr::arrange(chr_order, start, end) %>%
  dplyr::select(-chr_clean, -chr_numeric, -chr_order, -start, -end, -chr)

removed_orf_ids <- anti_join(orfs, filtered_table, by = c("gene_id", "gene_name", "Protein", "transcript_id", "orf_id")) %>%
  pull(orf_id)

# FOR TESTING PURPOSES REMOVE LATER
# Write unfiltered harmonised table
orfs_sorted <- orfs %>%
  tidyr::separate(ORF_ranges, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  tidyr::separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  dplyr::mutate(
      chr_clean = gsub("^chr", "", chr),
      chr_numeric = suppressWarnings(as.numeric(chr_clean)),
      chr_order = ifelse(is.na(chr_numeric), chr_clean, chr_numeric)
  ) %>%
  dplyr::arrange(chr_order, start, end) %>%
  dplyr::select(-chr_clean, -chr_numeric, -chr_order, -start, -end, -chr)
write.table(orfs_sorted, file = "unfiltered_harmonised_table.csv",
            sep = ",",
            quote = F,
            row.names = F)

# Write output to csv file
write.table(filtered_table_sorted, file = "harmonised_table.csv",
            sep = ",",
            quote = F,
            row.names = F)

# Save removed orf_ids to a text file
writeLines(removed_orf_ids, "removed_orf_ids.txt")

