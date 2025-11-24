#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(GenomicFeatures)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
orf_file <- args[1]
reference_gtf <- args[2]

# Import ORF gtf file
orf_gtf <- import(orf_file)

# Import reference GTF
ref_gtf <- rtracklayer::import(reference_gtf) %>% as.data.frame()

# Obtain reference exons
txdb <- txdbmaker::makeTxDbFromGFF(reference_gtf)
ref_exons <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Add 3 to CDS taking into account exon boundaries

# Group ORFs by ORF_id
orf_gtf <- orf_gtf[orf_gtf$type == "CDS"]

# Clean transcript IDs: original logic for ORFs, second-underscore for TCONS
mcols(orf_gtf)$transcript_id <- sapply(mcols(orf_gtf)$transcript_id, function(id) {
  if (grepl("^TCONS_", id) || grepl("^TM_", id)) {
    # keep everything before the second underscore
    parts <- unlist(strsplit(id, "_"))
    paste(parts[1:2], collapse = "_")
  } else {
    # keep everything before the first underscore
    sub("_.*$", "", id)
  }
})

# Group ORFs by ORF_id
orf_gtf_group <- split(orf_gtf, orf_gtf$ORF_id)

# Get transcript list
orf_transcripts <- data.frame(orf_gtf) %>% 
  distinct(ORF_id, transcript_id) %>% 
  arrange(ORF_id) %>% 
  pull(transcript_id)

# Map ORFs to transcript exons
tx_coord <- pmapToTranscripts(orf_gtf_group, ref_exons[orf_transcripts]) 
tx_coord <- unlist(tx_coord)
end(tx_coord) <- end(tx_coord) + 3 # Add plus 3 to end coordinate strand aware

# Map ORFs back to genome to obtain original coordinates
genomic_coord <- pmapFromTranscripts(tx_coord, ref_exons[orf_transcripts]) %>% 
  setNames(names(orf_gtf_group)) %>% 
  data.frame() %>% 
  dplyr::filter(hit) %>% 
  group_by(orf_id = group_name, chr = seqnames, strand) %>% 
  summarise(starts = paste0(start, collapse = ","),
            ends = paste0(end, collapse = ","),
            genomic_start = min(start),
            genomic_end = max(end),
            .groups = "keep")

# Split into one row for each start and stop pair
genomic_coord_expanded <- genomic_coord %>% 
  separate_rows(starts, ends, sep = ",") %>% 
  mutate(starts = as.numeric(starts),
         ends = as.numeric(ends)) %>%
  dplyr::rename(ORF_id = orf_id) %>%
  ungroup() %>%
  dplyr::select(ORF_id, starts, ends) 

# Obtain reference transcript features
tx_info <- ref_gtf %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, transcript_biotype, gene_id)

# Obtain reference gene features
gene_info <- ref_gtf %>%
  dplyr::filter(type == "gene") %>%
  dplyr::select(gene_id, gene_biotype, gene_name)

# Turn ORF rows into gtf ready GRanges object
orf_gtf_df <- as.data.frame(orf_gtf) %>%
  dplyr::filter(type == "CDS") %>%
  dplyr::select(seqnames, strand, source, type, transcript_id, ORF_id) %>%
  distinct() %>%
  # Join transcript biotype
  left_join(tx_info, by = "transcript_id") %>%
  # Join gene biotype
  left_join(gene_info, by = "gene_id") %>%
  # Join expanded genomic coordinates per ORF
  left_join(genomic_coord_expanded, by = "ORF_id") %>%
  mutate(
    start = as.integer(starts),
    end   = as.integer(ends)
  ) %>%
  dplyr::select(-starts, -ends) %>%
  makeGRangesFromDataFrame(
    seqnames.field = "seqnames",
    start.field    = "start",
    end.field      = "end",
    strand.field   = "strand",
    keep.extra.columns = TRUE
  )

# write proper GTF
rtracklayer::export(orf_gtf_df, "RiboTIE.gtf", format = "gtf")
