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
    length_kb = n() / 1000
  )

#intersect_files <- grep(
#  invert = T,
#  x = list.files(bedfile_loc, pattern = "intersect.bed", full.names = T),
#  value = T
#)

print(bedfile_loc)

# Split the string into a vector of filenames
bed_file_list <- strsplit(bedfile_loc, " ")[[1]]
print(bed_file_list)


#bedfile_loc <- getwd()  # Set to the current working directory
#print(bedfile_loc)
#intersect_files <- grep("intersect.bed", list.files(bedfile_loc), invert = FALSE, value = TRUE)
#print(intersect_files)

# Initiate DFs for populating
ppm <- data.frame(row.names = ref_ORFs_codons$ref_id)
psites <- data.frame(row.names = ref_ORFs_codons$ref_id)

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
    dplyr::summarize(psites = sum(score)) %>%
    dplyr::mutate(psitesfrac = psites / tot)
  
  # Select pX with highest fraction of p-sites
  p_select <- psites_test %>%
    dplyr::group_by(ref_id) %>%
    dplyr::filter(psitesfrac == max(psitesfrac)) %>%
    dplyr::transmute(frame_select = ref_frame)
  
  # For each ORF, only include p-sites from selected pX
  intersect_bed <- left_join(intersect_bed, p_select) %>%
    dplyr::group_by(ref_id) %>%
    dplyr::filter(ref_frame == frame_select)
  
  # Calculate p-sites
  psites_overlap <- intersect_bed %>%
    dplyr::group_by(ref_id) %>%
    dplyr::summarize(psites = sum(score)) %>%
    dplyr::full_join(ref_ORFs_codons) %>%
    dplyr::mutate(psites = ifelse(is.na(psites), 0, psites)) %>%
    dplyr::mutate(psites_perkb = psites / length_kb)
  
  # Calculate PPM
  scaling_factor <- sum(psites_overlap$psites_perkb) / 1000000
  psites_overlap$ppm <- psites_overlap$psites_perkb / scaling_factor
  
  rownames(psites_overlap) <- psites_overlap$ref_id
  psites_overlap <- psites_overlap[rownames(psites), ]
  
  ppm[[sample_name]] <- psites_overlap$ppm
  psites[[sample_name]] <- psites_overlap$psites
}

write.table(
  file = paste0(analysis_name, "_psites_permillion.txt"),
  x = ppm,
  row.names = T,
  quote = F
)
write.table(
  file = paste0(analysis_name, "_psites.txt"),
  x = psites,
  row.names = T,
  quote = F
)
