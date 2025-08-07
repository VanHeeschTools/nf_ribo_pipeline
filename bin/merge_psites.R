#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ORFquant)
  library(dplyr)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
psites_file <- args[1]

# Function to process and merge p-sites from GRanges objects
processAndMergePSites <- function(p_sites_list) {
  # Separate + and - stranded p-sites
  p_sites_plus <- p_sites_list[strand(p_sites_list) == "+"]
  p_sites_min <- p_sites_list[strand(p_sites_list) == "-"]

  cov_plus <- coverage(p_sites_plus, weight = S4Vectors::mcols(p_sites_plus)$score)
  cov_plus <- GenomicRanges::GRanges(cov_plus)
  cov_plus <- cov_plus[cov_plus$score > 0]

  cov_min <- coverage(p_sites_min, weight = S4Vectors::mcols(p_sites_min)$score)
  cov_min <- GenomicRanges::GRanges(cov_min)
  cov_min <- cov_min[cov_min$score > 0]

  strand(cov_plus) <- "+"
  strand(cov_min) <- "-"

  p_sites_merged <- sort(c(cov_plus, cov_min))
  rm(p_sites_plus, p_sites_min, cov_plus, cov_min)
  return(p_sites_merged)
}

# Function to merge junctions (summing reads and unique_reads)
mergeJunctions <- function(junctions_list) {
  summed_reads <- S4Vectors::mcols(junctions_list[[1]])$reads
  summed_unique_reads <- S4Vectors::mcols(junctions_list[[1]])$unique_reads
  
  for (i in 2:length(junctions_list)) {
    summed_reads <- summed_reads + S4Vectors::mcols(junctions_list[[i]])$reads
    summed_unique_reads <- summed_unique_reads + S4Vectors::mcols(junctions_list[[i]])$unique_reads
  }
  
  junctions_merged <- junctions_list[[1]]
  S4Vectors::mcols(junctions_merged)$reads <- summed_reads
  S4Vectors::mcols(junctions_merged)$unique_reads <- summed_unique_reads
  
  return(junctions_merged)
}

# Load for_ORFquant files
for_orfquant_files <- read.delim(psites_file, header = F)[, 1]
num_samples <- length(for_orfquant_files)

if (num_samples == 1) {
  message("Only one sample detected. Loading for_ORFquant object and saving without merging.")
  for_orfquant_file <- get(load(for_orfquant_files[1]))  # Load for_ORFquant object
  # Assuming the object is named for_ORFquant or similar in that file
  save(for_orfquant_file, file = "Merged_for_ORFquant")
  
} else if (num_samples > 1) {

  p_sites_all_list <- GenomicRanges::GRangesList()
  p_sites_uniq_list <- GenomicRanges::GRangesList()
  p_sites_uniq_mm_list <- GenomicRanges::GRangesList()
  junctions_list <- list()

  for (i in seq_along(for_orfquant_files)) {
    print(basename(for_orfquant_files[i]))
    for_orfquant_file <- get(load(for_orfquant_files[i])) # Assuming this loads for_orfquant_file
    p_sites_all_list[[i]] <- for_orfquant_file$P_sites_all
    p_sites_uniq_list[[i]] <- for_orfquant_file$P_sites_uniq
    p_sites_uniq_mm_list[[i]] <- for_orfquant_file$P_sites_uniq_mm
    junctions_list[[i]] <- for_orfquant_file$junctions
    # Clean up
    rm(for_orfquant_file)
    rm(for_ORFquant)
    gc()
  }

  # Unlist, process, and merge for each P sites type
  merged_p_sites_all <- processAndMergePSites(unlist(p_sites_all_list, recursive = FALSE))
  merged_p_sites_uniq <- processAndMergePSites(unlist(p_sites_uniq_list, recursive = FALSE))
  merged_p_sites_uniq_mm <- processAndMergePSites(unlist(p_sites_uniq_mm_list, recursive = FALSE))

  merged_junctions <- mergeJunctions(junctions_list)

  # Initiate and populate output list in the style ORFquant requires
  for_ORFquant <- list()
  for_ORFquant$P_sites_all <- merged_p_sites_all
  for_ORFquant$P_sites_uniq <- merged_p_sites_uniq
  for_ORFquant$P_sites_uniq_mm <- merged_p_sites_uniq_mm
  for_ORFquant$junctions <- merged_junctions

  save(for_ORFquant, file = "Merged_for_ORFquant")
} else {
  stop("No files found in merge_psites step.")
}
