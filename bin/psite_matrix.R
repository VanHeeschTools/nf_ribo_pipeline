#!/usr/bin/env Rscript

# Load libraries 
suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(rtracklayer)
  library(stringr)
})

# Obtain input arguments
args <- commandArgs(trailingOnly = TRUE)
ref_bed <- args[1]
bedfile_loc <- args[2]
analysis_name <- "orf_table"

# Define frame
frame_map <- c(p0 = 0L, p1 = 1L, p2 = 2L)

# Get lengths of reference ORFs
ref_bed <- data.table::fread(ref_bed,
  col.names = c("chrom", "start", "end", "ref_id", "frame", "strand", "nt_position"))

ref_ORFs_codons <- ref_bed %>%
  dplyr::group_by(ref_id) %>%
  dplyr::summarize(
    n_codons  = n(),
    length    = n() * 3,
    length_kb = length / 1000,
    strand    = dplyr::first(strand),
    .groups = "drop"
  )

# Split the string into a vector of filenames
bed_file_list <- strsplit(bedfile_loc, " ")[[1]]

# Sort bedfile vector numerically-aware
bed_file_list <- str_sort(bed_file_list, numeric = TRUE)

# Define dataframes
ppm <- data.frame(orf_id = ref_ORFs_codons$ref_id)
psites <- data.frame(orf_id = ref_ORFs_codons$ref_id)
psites_all_frames <- data.frame(orf_id = ref_ORFs_codons$ref_id)

# Define vector of all input bedfiles
codon_accum_list <- vector("list", length(bed_file_list))

# Loops over all files, extracts P-sites
for (i in seq_along(bed_file_list)) {
  int_file <- bed_file_list[i]
  sample_name <- gsub(pattern = "_intersect.bed",
                      replacement = "",
                      x = basename(int_file))

  intersect_bed <- data.table::fread(
    int_file,
    col.names = c(
      "chrom", "start", "end", "transcript_id", "score", "strand",
      "chrom_ref", "start_ref", "end_ref", "ref_id", "ref_frame",
      "ref_strand", "nt_position"
    )
  )

  # Go through intersect bedfiles individually
  codon_accum_list[[i]] <- intersect_bed[
    , .(psites = sum(score)),
    by = .(ref_id, codon = ceiling(nt_position / 3), frame = frame_map[ref_frame])
  ]

  # Calculate p-sites (frame 0 only)
  psites_overlap <- intersect_bed %>%
    dplyr::filter(ref_frame == "p0") %>%
    dplyr::group_by(ref_id) %>%
    dplyr::summarize(psites = sum(score), .groups = "drop") %>%
    dplyr::full_join(ref_ORFs_codons, by = "ref_id") %>%
    dplyr::mutate(psites = ifelse(is.na(psites), 0, psites)) %>%
    dplyr::mutate(psites_perkb = psites / length_kb)

  # Calculate p-sites across all frames (no ref_frame filter)
  psites_overlap_allframes <- intersect_bed %>%
    dplyr::group_by(ref_id) %>%
    dplyr::summarize(psites = sum(score), .groups = "drop") %>%
    dplyr::full_join(ref_ORFs_codons, by = "ref_id") %>%
    dplyr::mutate(psites = ifelse(is.na(psites), 0, psites)) %>%
    dplyr::mutate(psites_perkb = psites / length_kb)

  if (nrow(psites_overlap) > 0) {
    scaling_factor <- sum(psites_overlap$psites_perkb) / 1000000
    psites_overlap$ppm <- psites_overlap$psites_perkb / scaling_factor

    psites_overlap <- psites_overlap %>%
      dplyr::rename(orf_id = ref_id)

    ppm <- dplyr::left_join(ppm, psites_overlap %>% dplyr::select(orf_id, ppm), by = "orf_id")
    psites <- dplyr::left_join(psites, psites_overlap %>% dplyr::select(orf_id, psites), by = "orf_id")

    colnames(ppm)[ncol(ppm)] <- sample_name
    colnames(psites)[ncol(psites)] <- sample_name
  } else {
    warning(paste("No data available for sample:", sample_name))
  }

  if (nrow(psites_overlap_allframes) > 0) {
    psites_overlap_allframes <- psites_overlap_allframes %>%
      dplyr::rename(orf_id = ref_id)

    psites_all_frames <- dplyr::left_join(
      psites_all_frames,
      psites_overlap_allframes %>% dplyr::select(orf_id, psites),
      by = "orf_id"
    )

    colnames(psites_all_frames)[ncol(psites_all_frames)] <- sample_name
  } else {
    warning(paste("No data available (all frames) for sample:", sample_name))
  }
}

# Combine psites across all samples
codon_accum <- data.table::rbindlist(codon_accum_list)[
  , .(psites = sum(psites)), by = .(ref_id, codon, frame)
]
rm(codon_accum_list)

orf_len <- ref_ORFs_codons %>%
  dplyr::select(ref_id, orf_len = n_codons, strand)

codon_accum <- codon_accum %>%
  dplyr::left_join(orf_len, by = "ref_id")

# PIF: % of P-sites in frame 0, out of all P-sites in the ORF
pif <- codon_accum %>%
  dplyr::count(ref_id, frame, wt = psites, name = "psites") %>%
  tidyr::complete(ref_id, frame = 0:2, fill = list(psites = 0)) %>%
  tidyr::pivot_wider(names_from = frame, values_from = psites, names_prefix = "f") %>%
  dplyr::mutate(
    total_psites = f0 + f1 + f2,
    pif = dplyr::if_else(total_psites > 0, f0 / total_psites * 100, 0)
  )

# Uniformity: fraction of the ORF's codons where frame 0 is the majority (>1/3)
setDT(codon_accum)
codon_totals <- codon_accum[
  , .(f0 = sum(psites[frame == 0]), total = sum(psites)),
  by = .(ref_id, codon)
]
codon_totals[, pass := total > 0 & f0 / total > 1/3]
n_pass_dt <- codon_totals[, .(n_pass = sum(pass)), by = ref_id]

uniformity <- n_pass_dt %>%
  dplyr::left_join(orf_len, by = "ref_id") %>%
  dplyr::mutate(uniformity = dplyr::coalesce((n_pass / orf_len) * 100, 0)) %>%
  dplyr::select(ref_id, uniformity)

# Join statistics into final table
scores <- orf_len %>%
  dplyr::select(orf_id = ref_id, strand, orf_len_with_stop = orf_len) %>%
  dplyr::left_join(pif %>% dplyr::select(orf_id = ref_id, f0, f1, f2, total_psites, pif), by = "orf_id") %>%
  dplyr::left_join(uniformity %>% dplyr::rename(orf_id = ref_id), by = "orf_id") %>%
  dplyr::mutate(
    dplyr::across(c(f0, f1, f2, total_psites), ~ dplyr::coalesce(., 0)),
    pif        = round(dplyr::coalesce(pif, 0), 2),
    uniformity = round(dplyr::coalesce(uniformity, 0), 2)
  ) %>%
  dplyr::relocate(strand, orf_len_with_stop, .after = dplyr::last_col())

# Write outputs 
write.csv(file = paste0(analysis_name, "_psites_permillion.csv"), x = ppm, quote = FALSE, row.names = FALSE)
write.csv(file = paste0(analysis_name, "_psites_p0.csv"), x = psites, quote = FALSE, row.names = FALSE)
write.csv(file = paste0(analysis_name, "_psites_all_frames.csv"), x = psites_all_frames, quote = FALSE, row.names = FALSE)
write.csv(file = paste0(analysis_name, "_translation_scores.csv"), x = scores, quote = FALSE, row.names = FALSE)