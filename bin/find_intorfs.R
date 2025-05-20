#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(magrittr)
    library(dplyr)
    library(rtracklayer)
    library(GenomicRanges)
    library(data.table)
)}

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
transcriptome_gtf <- args[1]
transcriptome_p_sites <- args[2]
orfquant_annotation <- args[3]
psites_orfquant <- args[4]
price_annotation <- args[5]
psites_price <- args[6]

output_file <- "orf_id_canonical_overlaps.csv"


process_psites_transcriptome_gtf <- function(gtf_file, bed_file) {
  # Load and filter GTF annotation
  cds_annotation <- rtracklayer::import(gtf_file) %>%
    as.data.frame() %>%
    filter(type == "CDS") %>%
    dplyr::select(transcript_id, gene_id) %>%
    distinct()
  
  # Read BED file using fread()
  psites_cds <- fread(bed_file, header = FALSE, sep = "\t", select = c(1:6)) %>%
    as.data.frame()  # Convert to data.frame for compatibility
  
  # Rename columns
  colnames(psites_cds) <- c("seqnames", "start", "stop", "transcript_id", "codon_pos", "strand")
  
  # Filter rows to reduce dataset size
  psites_cds <- psites_cds %>%
    filter(codon_pos == "p0") %>%  
    semi_join(cds_annotation, by = "transcript_id")
  # Convert to GRanges object
  cds_grl <- psites_cds %>%
    left_join(cds_annotation, by = "transcript_id") %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  return(cds_grl)
}

process_psites_orfcaller <- function(orfcaller, annotation_file, bed_file) {
  # Load the annotation based on ORF caller
  if (orfcaller == "orfquant") {
    annotation <- fread(annotation_file) %>%
      as.data.frame() %>%
      dplyr::select(ORF_id_tr) %>%
      distinct()
    id_column <- "ORF_id_tr"  # Column to filter by
  } else if (orfcaller == "price") {
    annotation <- fread(annotation_file) %>%
      as.data.frame() %>%
      dplyr::select(name) %>%
      distinct()
    id_column <- "name"  # Column to filter by
  } else {
    stop("Invalid ORF caller")
  }
  
  # Read BED file
  psites <- fread(bed_file, header = FALSE, sep = "\t", select = c(1:6)) %>%
    as.data.frame()
  
  # Rename columns for clarity
  colnames(psites) <- c("seqnames", "start", "stop", "orf_id", "codon_pos", "strand")
  
  # Filter rows based on ORF caller
  psites_filtered <- psites %>%
    filter(codon_pos == "p0") %>%  # Filter before join to speed up function
    semi_join(annotation, by = c("orf_id" = id_column))

  return(psites_filtered)
}

cds_grl_test <- process_psites_transcriptome_gtf(transcriptome_gtf, transcriptome_p_sites)



# Obtain intorfs from orfcallers
orfquant_intorfs <- process_psites_orfcaller("orfquant", 
                                             orfquant_annotation, 
                                             psites_orfquant) 

price_intorfs <- process_psites_orfcaller("price", 
                                          price_annotation, 
                                          psites_price)

# Combine intorf tables
intorfs <- rbind(price_intorfs, orfquant_intorfs) %>%
  GenomicRanges::makeGRangesFromDataFrame()
names(intorfs) <- c(price_intorfs$orf_id, orfquant_intorfs$orf_id)

# Compare intORF p0 with all CDS P-sites
ol <- GenomicRanges::findOverlaps(query = cds_grl,
                                  subject = intorfs)

hit_orf_ids <- unique(names(intorfs[subjectHits(ol), ]))


write.csv(hit_orf_ids, output_file, row.names = FALSE)

