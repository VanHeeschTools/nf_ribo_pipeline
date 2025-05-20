#!/usr/bin/env Rscript

# This script fixes the ORFquant GTF output
# It has some incorrect ID values for the ORFs

suppressPackageStartupMessages({
    library(dplyr)
    library(ORFquant)
    library(Biostrings)
    library(tibble)
    library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
orfquant_results <- args[1]
rannot <- args[2]
orfquant_prefix <- args[3]
package_install_loc <- args[4]

paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Get workdir
workdir <- getwd()  

# Load files
orfquant_orfs <- get(load(orfquant_results))
load_annotation(rannot)

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
ORFs_readthroughs <- ORFquant_results$ORFs_readthroughs

# Fix protein FASTA
if (is(ORFs_readthroughs$Protein, "list")) {
  proteins_readthrough <- Biostrings::AAStringSet(lapply(ORFs_readthroughs$Protein, "[[", 1))
}else {
  proteins_readthrough <- Biostrings::AAStringSet(ORFs_readthroughs$Protein)
}
if (length(proteins_readthrough) > 0) {
  names(proteins_readthrough) <- paste(
    ORFs_readthroughs$ORF_id_tr,
    ORFs_readthroughs$gene_biotype,
    ORFs_readthroughs$gene_id,
    "readthrough",
    "readthrough",
    sep = "|"
  )
  proteins_readthrough <- GenomicRanges::narrow(proteins_readthrough, start = start(proteins_readthrough)[1] +
                                                  1)
  proteins_readthrough <- Biostrings::AAStringSet(gsub(
    proteins_readthrough,
    pattern = "[*]",
    replacement = "X"
  ))
}

proteins <- Biostrings::AAStringSet(ORFs_tx$Protein)
names(proteins) <- paste(
  ORFs_tx$ORF_id_tr,
  ORFs_tx$gene_biotype,
  ORFs_tx$gene_id,
  ORFs_tx$ORF_category_Gen,
  ORFs_tx$ORF_category_Tx_compatible,
  sep = "|"
)
proteins <- c(proteins, proteins_readthrough)
Biostrings::writeXStringSet(proteins, filepath = paste(orfquant_prefix, "Protein_sequences_fixed.fasta", sep = "_"))

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


# Export output files
filepath_fasta <- file.path(workdir, paste(orfquant_prefix, "Protein_sequences_fixed.fasta", sep = "_"))
Biostrings::writeXStringSet(proteins, filepath = filepath_fasta)


# EXTEND LAST CDS +3bp
cds <- all[all$type == "CDS"]

# Create tibble to track strand and positions
idxs <- tibble(
  idx = seq_along(cds),
  ORF = cds$ORF_id,
  s = start(cds),
  e = end(cds),
  strand = as.character(strand(cds))
) %>%
  group_by(ORF) %>%
  filter((strand == "+" & e == max(e)) | (strand == "-" & s == min(s))) %>%
  ungroup()

# Apply CDS extension
end(cds)[idxs$idx[idxs$strand == "+"]] <- end(cds)[idxs$idx[idxs$strand == "+"]] + 3L
start(cds)[idxs$idx[idxs$strand == "-"]] <- start(cds)[idxs$idx[idxs$strand == "-"]] - 3L

# Replace updated CDS in full GTF object
all[all$type == "CDS"] <- cds

# Export updated GTF
filepath_gtf <- file.path(workdir, paste(orfquant_prefix, "Detected_ORFs_fixed.gtf", sep = "_"))
suppressWarnings(rtracklayer::export.gff2(object = all, con = filepath_gtf))


