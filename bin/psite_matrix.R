#!/usr/bin/env Rscript

# Load libraries 
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(data.table)
  library(rtracklayer)
})

# Obtain input arguments
args = commandArgs(trailingOnly = TRUE)
ref_bed <- args[1]
bedfile_loc <- args[2]
analysis_name <- "orf_table"

# Get lengths of reference ORFs 
ref_bed <- data.table::fread(ref_bed,
    col.names = c("chrom", "start", "end", "ref_id", "frame", "strand", "nt_position")) 

ref_ORFs_codons <- ref_bed %>%
  dplyr::group_by(ref_id) %>%
  dplyr::summarize(
    n_codons = n(),
    length = n() * 3,
    length_kb = length / 1000,
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
      "ref_strand",
      "nt_position"
    )
  )
  
  # Calculate p-sites
  psites_overlap <- intersect_bed %>%
    dplyr::filter(ref_frame == "p0") %>% 
    dplyr::group_by(ref_id) %>%
    dplyr::summarize(psites = sum(score), .groups = "drop") %>%
    dplyr::full_join(ref_ORFs_codons, by = "ref_id") %>%
    dplyr::mutate(psites = ifelse(is.na(psites), 0, psites)) %>%
    dplyr::mutate(psites_perkb = psites / length_kb)
  
  # Ensure psites_overlap has the correct number of rows
  # TODO: Check that PPM are not counted double
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
