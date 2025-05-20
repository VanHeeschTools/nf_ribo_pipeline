#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(magrittr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
})

transcriptome_p_sites <-
psites_orfquant <-
psites_price <-
output_file <-


#' Obtain p0_sites from input bed file in order to find inframe ORFs
#' @param column_id, string: name value the id column should have
#' @param bed_file, input bed file with p0 sites 
#'
#' @return All p0 coords from CDS or ORFs, GRanges object with p0 site location
process_psites <- function(column_id, bed_file) {

  # Read BED file using fread()
  p0_sites <- data.table::fread(bed_file, sep = "\t", header = FALSE, select = 1:6,
                              col.names = c("seqnames", 
                                            "start", 
                                            "stop", 
                                            column_id, 
                                            "codon_pos", 
                                            "strand")) %>%
    dplyr::filter(codon_pos == "p0") # Only keep the p0 rows
  
  # Convert to GRanges object
  filtered_output <- p0_sites %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  return(filtered_output)
}

# Step_1: Call the process_psites function for the CDS
column_id <- "transcript_id"
cds_grl <- process_psites(column_id, transcriptome_p_sites)

# Step_2: Call the process_psites function for the ORFcallers
column_id <- "orf_id"
orfquant_intorfs <- process_psites(column_id, psites_orfquant) 
price_intorfs <- process_psites(column_id, psites_price)

# Step_3: Combine intORF tables
# Ignore the sequence level warning
# Not all contigs have an ORF so sequences might be missing
intorfs <- suppressWarnings(c(price_intorfs, orfquant_intorfs))  
names(intorfs) <- c(as.character(price_intorfs$orf_id),   
                    as.character(orfquant_intorfs$orf_id))

# Step_4: Compare intORF p0 with all CDS P-sites
ol <- GenomicRanges::findOverlaps(query = cds_grl,
                                  subject = intorfs)

# Step_5: Obtain all unique ORF ids that are in the same frame as a CDS
hit_orf_ids <- unique(names(intorfs[subjectHits(ol), ]))