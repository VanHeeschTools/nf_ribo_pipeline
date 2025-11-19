#!/usr/bin/env Rscript

# This script fixes the ORFquant GTF output
# It has some incorrect ID values for the ORFs
# It also add 3 bp to the last coords to take the stop codon into account

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ORFquant)
  library(Biostrings)
  library(tibble)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
})

# Obtain input data from arguments
args <- commandArgs(trailingOnly = TRUE)
orfquant_results <- args[1]
rannot <- args[2]
package_install_loc <- args[3]
reference_gtf <- args[4]

paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Get workdir
workdir <- getwd()  

# Load files
orfquant_orfs <- get(load(orfquant_results))
load_annotation(rannot)

# Import reference GTF
ref_gtf <- rtracklayer::import(reference_gtf) %>% as.data.frame()

# Obtain reference exons
txdb <- txdbmaker::makeTxDbFromGFF(reference_gtf)
ref_exons <- exonsBy(txdb, by = "tx", use.names = TRUE)

# Fix ORFquant annotation
ORFs_txs_feats <- orfquant_orfs$ORFs_txs_feats
selected_txs <- sort(unique(unlist(ORFs_txs_feats$txs_selected)))

ORFs_tx <- orfquant_orfs$ORFs_tx

# USE ORFS_tx to get the correct info per ORF
map_tx_genes <- S4Vectors::mcols(ORFs_tx)[, c(
  "ORF_id_tr",
  "gene_id",
  "gene_biotype",
  "gene_name",
  "transcript_id",
  "transcript_biotype",
  "P_sites",
  "ORF_pct_P_sites",
  "ORF_pct_P_sites_pN",
  "ORFs_pM"
)]

# Fix ORFs_gen based on ORFs_tx
ORFs_gen <- orfquant_orfs$ORFs_gen

match_ORF <- match(names(ORFs_gen), map_tx_genes$ORF_id_tr)
ORFs_gen$transcript_id <- map_tx_genes[match_ORF, "transcript_id"]
match_tx <- match(ORFs_gen$transcript_id, map_tx_genes$transcript_id)

ORFs_gen$transcript_id <- map_tx_genes[match_ORF, "transcript_id"]
ORFs_gen$transcript_biotype <- map_tx_genes[match_tx, "transcript_biotype"]
ORFs_gen$gene_id <- map_tx_genes[match_tx, "gene_id"]
ORFs_gen$gene_biotype <- map_tx_genes[match_tx, "gene_biotype"]
ORFs_gen$gene_name <- map_tx_genes[match_tx, "gene_name"]
ORFs_gen$ORF_id <- map_tx_genes[match_ORF, "ORF_id_tr"]
ORFs_gen$P_sites <- round(map_tx_genes[match_ORF, "P_sites"], digits = 4)
ORFs_gen$ORF_pct_P_sites <- round(map_tx_genes[match_ORF, "ORF_pct_P_sites"],
                                  digits = 4)
ORFs_gen$ORF_pct_P_sites_pN <- round(map_tx_genes[match_ORF, "ORF_pct_P_sites_pN"],
                                    digits = 4)
ORFs_gen$ORFs_pM <- round(map_tx_genes[match_ORF, "ORFs_pM"], 
                          digits = 4)
#ORFs_readthroughs <- ORFquant_results$ORFs_readthroughs


# Create new fixed ORFquant GTF
map_tx_genes <- GTF_annotation$trann
ORFs_gen$type = "CDS"
exs_gtf <- unlist(GTF_annotation$exons_txs[selected_txs])
S4Vectors::mcols(exs_gtf) <- NULL
exs_gtf$transcript_id <- names(exs_gtf)
exs_gtf$transcript_biotype <- map_tx_genes[match(exs_gtf$transcript_id, map_tx_genes$transcript_id), "transcript_biotype"]
exs_gtf$gene_id <- map_tx_genes[match(exs_gtf$transcript_id, map_tx_genes$transcript_id), "gene_id"]
exs_gtf$gene_biotype <- map_tx_genes[match(exs_gtf$transcript_id, map_tx_genes$transcript_id), "gene_biotype"]
exs_gtf$gene_name <- map_tx_genes[match(exs_gtf$transcript_id, map_tx_genes$transcript_id), "gene_name"]
S4Vectors::mcols(exs_gtf)[, names(S4Vectors::mcols(ORFs_gen))[!names(S4Vectors::mcols(ORFs_gen)) %in% names(S4Vectors::mcols(exs_gtf))]] <-
  NA
exs_gtf$type <- "exon"
S4Vectors::mcols(ORFs_gen) <- S4Vectors::mcols(ORFs_gen)[, names(S4Vectors::mcols(exs_gtf))]
all <- sort(c(exs_gtf, ORFs_gen))
all$`source` = "ORFquant"
names(all) <- NULL


# Add 3 to CDS taking into account exon boundaries

# Group ORFs by ORF_id
orf_gtf <- all[all$type == "CDS"]

# Clean ORF_ids to keep transcript_id
mcols(orf_gtf)$transcript_id <- sapply(mcols(orf_gtf)$transcript_id, function(id) {
  if (grepl("^TCONS_", id)) {
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
rtracklayer::export(orf_gtf_df, "ORFquant.gtf", format = "gtf")
