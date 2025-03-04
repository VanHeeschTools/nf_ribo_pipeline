#!/usr/bin/env Rscript

# Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(data.table)
  library(rtracklayer)
})

# Input args --------------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)

ref_bedf <- args[1]
bedfile_loc <- args[2]
analysis_name <- args[3]

# Get lengths of reference ORFs -------------------------------------------
ref_bed <- data.table::fread(ref_bedf,
                             col.names = c("chrom", "start", "end", "ref_id", "frame", "strand")) %>%
  subset(!grepl("pATG|pST", .$frame))

ref_ORFs_codons <- ref_bed %>%
  dplyr::group_by(ref_id) %>%
  dplyr::summarize(
    n_codons = n() / 3,
    length = n(),
    length_kb = n() / 1000,
    .groups = "drop"
  )

# Split the string into a vector of filenames
bed_file_list <- strsplit(bedfile_loc, " ")[[1]]

# Initiate DFs for populating (with ref_id as a column)
ppm <- data.frame(orf_id = ref_ORFs_codons$ref_id)
psites <- data.frame(orf_id = ref_ORFs_codons$ref_id)

# Loops over all files, extracts P-sites
for (int_file in bed_file_list) {
  sample_name = gsub(pattern = "_intersect.bed",
                     replacement = "",
                     x = basename(int_file))
  
  intersect_bed <- data.table::fread(
    int_file,
    col.names = c(
      "chrom",
      "start",
      "end",
      "transcript_id",
      "score",
      "strand",
      "chrom_ref",
      "start_ref",
      "end_ref",
      "ref_id",
      "ref_frame",
      "ref_strand"
    )
  )
  
  # Calculate fraction of p0, p1 and p2
  psites_test <- intersect_bed %>%
    # Skip start and stop codons
    subset(!grepl("pATG|pST", .$ref_frame)) %>%
    # Do the following calculations based on ref_id
    dplyr::group_by(ref_id) %>%
    dplyr::mutate(tot = sum(score)) %>%
    dplyr::group_by(ref_id, ref_frame, tot) %>%
    dplyr::summarize(psites = sum(score), .groups = "drop") %>%
    dplyr::mutate(psitesfrac = psites / tot)
  
  # Select pX with highest fraction of p-sites
  p_select <- psites_test %>%
    dplyr::group_by(ref_id) %>%
    dplyr::filter(psitesfrac == max(psitesfrac)) %>%
    dplyr::slice(1) %>%  # Ensure only one row per ref_id
    dplyr::transmute(frame_select = ref_frame)
  
  # For each ORF, only include p-sites from selected pX
  intersect_bed <- left_join(intersect_bed, p_select, by = "ref_id") %>%
    dplyr::group_by(ref_id) %>%
    dplyr::filter(ref_frame == frame_select)
  
  # Calculate p-sites
  psites_overlap <- intersect_bed %>%
    dplyr::group_by(ref_id) %>%
    dplyr::summarize(psites = sum(score), .groups = "drop") %>%
    dplyr::full_join(ref_ORFs_codons, by = "ref_id") %>%
    dplyr::mutate(psites = ifelse(is.na(psites), 0, psites)) %>%
    dplyr::mutate(psites_perkb = psites / length_kb)
  
  # Ensure psites_overlap has the correct number of rows
  if (nrow(psites_overlap) > 0) {
    # Calculate PPM
    scaling_factor <- sum(psites_overlap$psites_perkb) / 1000000
    psites_overlap$ppm <- psites_overlap$psites_perkb / scaling_factor
    
    # Keep orf_id as a column instead of rownames
    psites_overlap <- psites_overlap %>%
      dplyr::rename(orf_id = ref_id)
    
    # Safely join the new sample data to ppm and psites
    ppm <- dplyr::left_join(ppm, psites_overlap %>% dplyr::select(orf_id, ppm), by = "orf_id")
    psites <- dplyr::left_join(psites, psites_overlap %>% dplyr::select(orf_id, psites), by = "orf_id")
    
    # Rename columns to reflect sample names
    colnames(ppm)[ncol(ppm)] <- sample_name
    colnames(psites)[ncol(psites)] <- sample_name
  } else {
    warning(paste("No data available for sample:", sample_name))
  }
}

# Write outputs to CSV, with orf_id as a column
write.csv(
  file = paste0(analysis_name, "_psites_permillion.csv"),
  x = ppm,
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  file = paste0(analysis_name, "_psites.csv"),
  x = psites,
  quote = FALSE,
  row.names = FALSE
)
