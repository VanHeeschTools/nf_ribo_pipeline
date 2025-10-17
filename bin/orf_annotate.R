#!/usr/bin/env Rscript

# =============================================================================
# 01 | LOAD LIBRARIES ----
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(AnnotationDbi)
  library(Biostrings)
  library(parallel)
  library(GenomicFeatures)
  library(S4Vectors)
  library(txdbmaker)
})

# =============================================================================
# 02 | LOAD DATA ----
#   * Read input files
#   * Load custom BSgenome package
#   * Create variables required for analysis
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
orfs_loc <- args[1]
orfcaller_bed_file <- args [2]
ref_bed_file <- args[3]
gtf <- args[4]
orfcaller <- args[5]
package_install_loc <- args[6]
bsgenome_path <- args[7]

# Load custom BSgenome package
paths <- c(package_install_loc, .libPaths())
.libPaths(paths)
bsgenome_path = bsgenome_path
bsgenome_name <- basename(bsgenome_path)
#bsgenome_dirs <- list.dirs(bsgenome_path, full.names = FALSE, recursive = FALSE)
#bsgenome_name <- bsgenome_dirs[grepl("^BSgenome\\.", bsgenome_dirs)]
load_bsgenome <- function(pkg) {
  ns <- asNamespace(pkg)
  bsgenome_obj <- Filter(function(x) inherits(get(x, envir = ns), "BSgenome"), ls(ns))
  get(bsgenome_obj[[1]], envir = ns)
}
genome <- load_bsgenome(bsgenome_name)


# Extract CDS regions of annotated genes from txdb
txdb <- txdbmaker::makeTxDbFromGFF(gtf, format = "gtf")
transcript_ids <- AnnotationDbi::keys(txdb, keytype = "TXNAME")

# Obtain all genes and their corresponding CDS coords
cds_gene <- GenomicFeatures::cdsBy(txdb, "gene")
cds_gene_unlist <- unlist(cds_gene)

# Obtain all transcript ids and their corresponding gene ids
tx2gene <- AnnotationDbi::select(txdb, transcript_ids, "GENEID", "TXNAME")

# Import the GTF file as data frame
gtf_df <- as.data.frame(rtracklayer::import(gtf)) %>%
  filter(type == "gene") %>%
  dplyr::select(GENEID = gene_id, 
                GENEBIOTYPE = gene_biotype,
                GENENAME = gene_name)  

# =============================================================================
# 03 | FUNCTIONS ----
#   * Define all functions
# =============================================================================

#' Prepare PRICE ORF Data for Analysis
#'
#' Imports and processes ORFs predicted by PRICE from a BED file, generating
#' both genomic ranges and ORF-level metadata.
#'
#' @param price_orfs_loc A character string indicating the file path to the 
#'   PRICE ORFs in BED format.
#' @param tx2gene A data frame mapping transcript IDs to gene IDs
#' @param gtf_df A data frame containing gene annotations from a GTF file
#'
#' @return A list with two elements:
#' orf_ranges A GRangesList object representing the genomic 
#'   coordinates (blocks) of the predicted ORFs.
#' price_orf_df A data frame containing metadata for each ORF, including 
#'   ORF ID, transcript and gene identifiers, strand, range, gene name
#'   and biotype, the predicted protein sequence and length are added.
#' 
prepare_price <- function(price_orfs_loc, tx2gene, gtf_df) {
  # This functions imports the PRICE ORFs into R
  # Returns a list with the ORF genomic ranges (orf_ranges)
  # and the metadata info per ORF (price_orf_df)
  
  # Import PRICE ORF definitions
  price_orfs <- rtracklayer::import.bed(price_orfs_loc)
  # Use BED file blocks to establish exons 
  orf_ranges <- rtracklayer::blocks(price_orfs) 
  
  # Fix gene IDs for transcript IDs (which are correct)
  price_orf_df <- as.data.frame(price_orfs) %>%
    dplyr::select(seqnames,start,end,width,strand,name) %>%
    dplyr::mutate(
                  ORF_ranges = paste0(seqnames, ":", start, "-", end),
                  gene_id = stringr::str_split_i(name, "__", i = 1),
                  transcript_id = stringr::str_split_i(name, "__", i = 2),
                  start_codon = stringr::str_split_i(name, "__", i = 5)) %>%
    dplyr::left_join(tx2gene, by = c("transcript_id" = "TXNAME")) %>%
    dplyr::select(!(gene_id)) %>%
    # Add gene id and biotype
    dplyr::left_join(gtf_df, by = "GENEID")%>%  
    dplyr::rename(gene_id = GENEID, 
                  gene_biotype = GENEBIOTYPE, 
                  gene_name = GENENAME, 
                  orf_id = name) 
 
  # Obtain protein sequence for PRICE ORFs
  # Split the ORF table by strand for simplicity
  price_plus <- price_orf_df[which(price_orf_df$strand == "+"),]$orf_id
  price_minus <- price_orf_df[which(price_orf_df$strand == "-"),]$orf_id
  
  # Sort the minus strand (if necessary)
  orf_minus_sorted <- GenomicRanges::sort(orf_ranges[which(names(orf_ranges) %in% price_minus),], decreasing = T)  
  orf_plus_names <- orf_ranges[which(names(orf_ranges) %in% price_plus),]
  
  # Define a function to translate sequences and return protein sequences
  get_protein_seq <- function(orf_names, genome) {
    lapply(Biostrings::getSeq(x = genome, names = orf_names), 
           function(seq) {
             Biostrings::toString(Biostrings::translate(Biostrings::DNAString(paste(seq, collapse = "")), genetic.code = GENETIC_CODE_PRICE))
           })
  }
  
  # Check if genetic code is needed here and what it does
  # Create a custom genetic code (if required)
  GENETIC_CODE_PRICE <- Biostrings::GENETIC_CODE
  attr(GENETIC_CODE_PRICE, "alt_init_codons") <- unique(price_orf_df$start_codon)
  
  # Get protein sequences for both strands
  price_protein_plus <- get_protein_seq(orf_plus_names, genome)
  price_protein_minus <- get_protein_seq(orf_minus_sorted, genome)
  Protein = c(unlist(price_protein_plus), unlist(price_protein_minus))
  Protein = gsub("\\*$", "", Protein)  # Remove trailing '*' from Protein column

  # Combine results into a data frame
  price_proteins <- data.frame(
    orf_id = c(names(price_protein_plus), names(price_protein_minus)),
    Protein = Protein,
    Protein_Length = nchar(Protein)   
  )
  
  # Join with the original dataframe
  price_orf_df <- price_orf_df %>% left_join(price_proteins, by = "orf_id")
  
  return(list(orf_ranges, price_orf_df))
}

#' Prepare ORFquant Data for Analysis
#'
#' Imports the ORFquant result object and processes its contents to extract 
#' ORF genomic ranges and associated ORF metadata in a structured format. 
#' Returns both the ORF ranges and a metadata table with transcript- and 
#' gene-level annotations.
#'
#' @param orfquant_orfs_loc A character string specifying the path to the 
#' ORFquant `.RData` file that contains the ORFquant object.
#'
#' @return A list with two elements:
#' @return orf_ranges A named list of GRanges objects, one per ORF, representing 
#'   the genomic coordinates of each ORF including adjusted stop codon positions.
#' @return orf_table A data frame with metadata for each ORF, including ORF ID, 
#'   genomic range string, strand, protein sequence and length, 
#'   gene and transcript identifiers, biotypes, ORF categories, 
#'   and P-site information.
#'
#' @details
#' The function adjusts the ORF ranges to include the stop codon (which is 
#' omitted in ORFquant output), computes protein lengths, and ensures all 
#' entries are distinct. It joins the genomic and transcript-level 
#' information using the ORF ID.
#'
prepare_orfquant <- function(orfquant_orfs_loc) {

  # Imports the ORFquant object and outputs a list of ORF ranges (orf_ranges), 
  # and the ORF metadata (orfs_table) linked to the ORF ranges
  
  # Load ORFquant file and select ORF definitions
  orfquant_orfs <- get(load(orfquant_orfs_loc))
  orfs_tx_df <- data.frame(
    S4Vectors::mcols(orfquant_orfs$ORFs_tx)[, c("ORF_id_tr", 
                                                "Protein", 
                                                "gene_id", 
                                                "gene_biotype", 
                                                "gene_name", 
                                                "transcript_id", 
                                                "transcript_biotype", 
                                                "ORF_category_Tx", 
                                                "ORF_category_Gen", 
                                                "P_sites_raw", 
                                                "P_sites_raw_uniq")]) %>%
  dplyr::mutate(
    # Add protein length column
    Protein_Length = nchar(as.character(Protein))  
  ) %>%
    dplyr::relocate(Protein_Length, .after = Protein) %>%
    dplyr::distinct()
  
  # Generate table of genomic ORF locations
  orf_table <- data.frame(orfquant_orfs$ORFs_gen) %>%
    dplyr::mutate(
      ORF_id_tr = names(orfquant_orfs$ORFs_gen),
      # Extract strand before grouping
      strand = as.character(strand(orfquant_orfs$ORFs_gen))  
    ) %>%
    dplyr::group_by(ORF_id_tr) %>%
    dplyr::mutate(
      # Add the stop codon to range as this is missing in ORFquant
      ORF_ranges = dplyr::case_when(
        
        strand == "+" ~ paste0(seqnames, ":", min(start), "-", max(end) + 3),
        strand == "-" ~ paste0(seqnames, ":", min(start) - 3, "-", max(end)),
        TRUE ~ paste0(seqnames, ":", min(start), "-", max(end))  
      )      
    ) %>%
    
    dplyr::select(c("ORF_id_tr", "ORF_ranges", "strand")) %>%
    dplyr::distinct() %>%
    dplyr::left_join(orfs_tx_df, by = c("ORF_id_tr")) %>% 
    # Rename column to orf_id
    dplyr::rename("orf_id" = "ORF_id_tr") 
  
  orf_ranges <- split(orfquant_orfs$ORFs_gen, names(orfquant_orfs$ORFs_gen))
  
  return(list(orf_ranges,
              orf_table))
}

prepare_ribotie <- function(ribotie_orfs_loc) {
  # Load RiboTIE merged ORF table
  orf_summary_tbl <- read.csv(ribotie_orfs_loc) %>%
    dplyr::rename("orf_id" = "ORF_id",
                  "Protein" = "protein_seq") %>%
    dplyr::mutate(
      start_coord = ifelse(strand == "+", TIS_coord, LTS_coord - 3 ),
      end_coord = ifelse(strand == "+", LTS_coord + 3, TIS_coord)
    ) %>%
    dplyr::mutate(ORF_ranges = paste0(seqname, ":", start_coord, "-", end_coord),
                  Protein_Length = nchar(Protein))

  # Parse CDS_coords into flat GRanges 
  cds_df <- orf_summary_tbl %>%
    dplyr::select(orf_id, CDS_coords) %>%
    dplyr::filter(!is.na(CDS_coords)) %>%
    dplyr::mutate(CDS_coords = str_split(CDS_coords, ";\\s*")) %>%
    tidyr::unnest(CDS_coords) %>%
    dplyr::mutate(
      seqname = str_extract(CDS_coords, "^[^:]+"),
      start   = as.integer(str_extract(CDS_coords, "(?<=:)[0-9]+")),
      end     = as.integer(str_extract(CDS_coords, "(?<=-)[0-9]+")),
      strand  = str_extract(CDS_coords, "(?<=\\()[+-](?=\\))")
    )
  
  # Construct GRanges object
  cds_gr <- GRanges(
    seqnames = cds_df$seqname,
    ranges   = IRanges(cds_df$start, cds_df$end),
    strand   = cds_df$strand
  )
  names(cds_gr) <- cds_df$orf_id
  
  # Split GRanges by ORF_id
  orf_ranges <- split(cds_gr, cds_df$orf_id)
  
  # Remove columns that are no longer required
  #orf_summary_tbl <- orf_summary_tbl %>%
  #  dplyr::select(-TIS_coord,-LTS_coord,-start_coord,-end_coord)
  
  # Only keep required columns
  orf_summary_tbl <- orf_summary_tbl %>%
  dplyr::select(orf_id,seqname,ORF_len,
                transcript_id,start_codon,
                strand,ORF_type,
                Protein,gene_id,
                gene_name, gene_biotype,
                ORF_ranges, Protein_Length,
                CDS_coords
  )

  # Return list
  return(list(
    orf_ranges = orf_ranges,
    orf_table = orf_summary_tbl
  ))
}

#' Match Sequence Level Style Between Genomic Ranges
#'
#' Ensures that the sequence naming style (e.g., "chr1" vs "1") of the ORF 
#' genomic ranges matches that of the reference annotation object.
#'
#' @param orf_ranges An object containing ORF genomic ranges whose sequence 
#'   style needs to be updated.
#' @param annotated_gen The reference object with the desired seqlevelsStyle
#'
#' @return A modified version of the input GRanges object where the style 
#'   matches that of the annotated_gen.
#'
check_annot_style <- function(orf_ranges, annotated_gen) {
  # Match seqlevels style of orf_ranges to that of annotated_gen
  GenomeInfoDb::seqlevelsStyle(orf_ranges) <- GenomeInfoDb::seqlevelsStyle(annotated_gen)[1]
  return(orf_ranges)
}

#' Match predicted ORFs to annotated CDS regions, choosing the best overlap or 
#' in the case of no overlap falling back on gene-level CDSs
#'
#' For each ORF in `orf_ranges`, this function identifies the annotated CDS 
#' with which it has the greatest genomic overlap. If an ORF has no direct 
#' overlap with any CDS, it is assigned all CDS regions of its parent gene 
#' Returns a `GRangesList` aligned to `orf_ranges`.
#'
#' @param orf_ranges A `GRanges` or `GRangesList` of predicted ORFs.
#' @param orf_table A data.frame with columns `orf_id` and `gene_id`.
#' @param annotated_gen A `GRanges` of annotated CDS regions.
#' @param annotated_gen_unlist A flattened `GRanges` of all CDS regions.
#'
#' @return A `GRangesList` with one element per ORF:
#'   - If the ORF overlaps any annotated CDS, that single best-overlapping CDS 
#'     (`GRanges`) is returned.  
#'   - If no overlaps are found, all CDSs of the ORF’s parent gene are returned.  
#'   - If the parent gene has no CDS, the element remains an empty `GRanges`.
check_orf_cds_similarity <- function(orf_ranges, orf_table, 
                                    annotated_gen, annotated_gen_unlist) {
  
  # Find all pairwise overlaps between predicted ORFs and annotated CDSs
  overlaps <- GenomicRanges::findOverlaps(orf_ranges, annotated_gen)
  
  # Compute the width (in base pairs) of each overlap region
  overlap_width <- sum(
    width(
      GenomicRanges::intersect(
        orf_ranges[S4Vectors::queryHits(overlaps)], 
        annotated_gen[S4Vectors::subjectHits(overlaps)]
      )
    )
  )

  # Build a table of ORF indices, CDS indices, and their overlap widths
  overlap_df <- data.frame(
    queryIdx     = S4Vectors::queryHits(overlaps),   # ORF index in orf_ranges
    subjectIdx   = S4Vectors::subjectHits(overlaps), # CDS index in annotated_gen
    overlapWidth = overlap_width,                    # measured overlap
    stringsAsFactors = FALSE
  )
  
  # Keep only the highest-overlap CDS for each ORF
  max_overlaps <- overlap_df[order(overlap_df$queryIdx, 
                                  -overlap_df$overlapWidth), ]
  max_overlaps <- max_overlaps[!duplicated(max_overlaps$queryIdx), ]
  
  # Extract ORF and CDS indices for these best overlaps
  query_idx   <- max_overlaps$queryIdx
  subject_idx <- max_overlaps$subjectIdx
  
  # Initialize a mapping of each ORF to its selected CDS (NA means no overlap)
  selected_overlaps <- data.frame(
    queryHits   = seq_along(orf_ranges),
    subjectHits = rep(NA_integer_, length(orf_ranges))
  )
  selected_overlaps$subjectHits[
    selected_overlaps$queryHits %in% query_idx
  ] <- subject_idx
  
  # Create an empty GRangesList to hold the result, one slot per ORF
  result_list <- GenomicRanges::GRangesList(
    rep(list(GenomicRanges::GRanges()), length(orf_ranges))
  )
  names(result_list) <- names(orf_ranges)
  
  # Fill in each ORF slot with its best-overlapping CDS, if any
  non_na_indices <- !is.na(selected_overlaps$subjectHits)
  result_list[selected_overlaps$queryHits[non_na_indices]] <-
    annotated_gen[selected_overlaps$subjectHits[non_na_indices]]
  
  # Identify ORFs that still have no assigned CDS (empty GRanges)
  no_overlap_idx   <- lengths(result_list) == 0
  no_overlap_names <- names(which(no_overlap_idx))
  
  # For ORFs with no CDS overlap, assign all CDSs of their parent gene
  result_list[no_overlap_idx] <- GenomicRanges::GRangesList(
    lapply(no_overlap_names, function(name) {
      # Lookup the parent gene for this ORF in the orf_table
      orf_parent_gene <- orf_table$gene_id[match(name, orf_table$orf_id)]
      # Return all CDS ranges from that gene
      annotated_gen_unlist[names(annotated_gen_unlist) == orf_parent_gene]
    })
  )
  
  return(result_list)
}

#' Compare CDS vs ORFcaller ORF, and make a new annotation based on a virtual 
#' longest possible ORF of the annotated set.  
#' @param orf_ranges Granges of orf_caller ORFs
#' @param orf_table metadata of orf_caller ORFs
#' @param cds_matches_grl ()
#' @param orf_caller String containing ORFCaller name. Required because 
#'.                  ORFquant does not include the STOP codon in the sequence
#' 
#' @return orf_table with new columns showing start codon, orf coord similarity 
#'         and the new orf annotation type.

# ORFquant does not include stops in the ORF, so we trim the CDS matches
annotate_new_orfs <- function(orf_ranges, orf_table, cds_matches_grl, orf_caller) {
  # Convert ORF caller name to lowercase for consistency
  orf_caller <- tolower(orf_caller)
  # Adjust stop codon handling based on ORF caller
  adjust_stop <- ifelse(orf_caller == "price", 0, 3)  
  
  # Calculate similarity between ORF and CDS range
  cds_range_similarity <- width(range(orf_ranges)) / (width(range(cds_matches_grl)) - adjust_stop)
  cds_strand <- ifelse(elementNROWS(cds_matches_grl) > 0, as.character(unique(strand(cds_matches_grl))), NA)
  
  # Extract ORF strand and determine start/stop positions
  orf_strand <- as.character(unlist(runValue(strand(orf_ranges))))
  orf_start <- ifelse(orf_strand == "+", min(start(orf_ranges)), max(end(orf_ranges)))
  orf_stop  <- ifelse(orf_strand == "+", max(end(orf_ranges)), min(start(orf_ranges)))
  
  # Determine annotated start and stop positions
  ann_start <- ifelse(cds_strand == "+", min(start(cds_matches_grl)), max(end(cds_matches_grl)))
  ann_stop  <- ifelse(cds_strand == "+", max(end(cds_matches_grl)) - adjust_stop, min(start(cds_matches_grl)) + adjust_stop)

  # Default ORF category as "Unknown"
  orf_category <- rep("Unknown", length(orf_ranges))
  
  # Function to classify ORFs based on start and stop positions
  pos_strand_idx <- orf_strand == "+"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start == ann_start] <- "ORF_annotated"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start < ann_start] <- "N_extension"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start > ann_start] <- "N_truncation"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop < ann_stop] <- "overl_uORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop < ann_start] <- "uORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop > ann_stop] <- "NC_extension"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop > ann_stop] <- "overl_dORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_stop & orf_stop > ann_stop] <- "dORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop < ann_stop] <- "intORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop < ann_stop] <- "C_truncation"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop > ann_stop] <- "C_extension"
  
  # Reverse start / stop for negative strand
  neg_strand_idx <- orf_strand == "-"
  orf_category[neg_strand_idx & orf_stop == ann_stop & orf_start == ann_start] <- "ORF_annotated"
  orf_category[neg_strand_idx & orf_stop == ann_stop & orf_start > ann_start] <- "N_extension"
  orf_category[neg_strand_idx & orf_stop == ann_stop & orf_start < ann_start] <- "N_truncation"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop > ann_stop] <- "overl_uORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop > ann_start] <- "uORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop < ann_stop] <- "NC_extension"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop < ann_stop] <- "overl_dORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start < ann_stop & orf_stop < ann_stop] <- "dORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop > ann_stop] <- "intORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop > ann_stop] <- "C_truncation"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop < ann_stop] <- "C_extension"
  
  # Mark novel ORFs if they do not match any CDS
  orf_category[lengths(cds_matches_grl) == 0] <- "novel"
  cds_range_similarity[lengths(cds_matches_grl) == 0] <- NA

  # Set start codon to ATG for ORFquant
  start_codon <- rep("ATG", length(orf_ranges))  # Default to ATG

  # Add start_codon to new_category_df
  new_category_df <- data.frame(
    orf_id = names(orf_ranges),
    cds_range_similarity = as.numeric(cds_range_similarity),
    orf_category_new = orf_category,
    start_codon = start_codon  # Use fixed start codon extraction
  ) 
  
  # Remove start_codon column for "price" ORF caller 
  # start_codon column already exists in price table
  if (tolower(orfcaller) == "orfquant"){
    orf_table <- orf_table %>%
      dplyr::left_join(new_category_df, by = "orf_id")
  }
  
  if (tolower(orfcaller) %in% c("price","ribotie")) {
    new_category_df <- new_category_df %>%
      dplyr::select(!start_codon)
    orf_table <- orf_table %>%
      dplyr::left_join(new_category_df, by = "orf_id")
  }
  
  return(orf_table)
}

#' Reclassify and Refine ORF Annotations Based on Gene Biotypes
#'
#' This function updates gene biotype categories into better fitting
#' groups and refines ORF category annotations (`orf_category_new`) based on
#' the updated gene biotypes and CDS similarity information.
#'
#' @param orf_table A data frame containing ORF annotations,
#'   including the columns `gene_biotype`, `orf_category_new`, and optionally
#'   `cds_range_similarity`.
#'
#' @return The input `orf_table` with two columns updated:
#'   - `gene_biotype`: grouped into better fitting categories.
#'   - `orf_category_new`: re-annotated based on the new `gene_biotype` and
#'      overlap with CDS (`cds_range_similarity`).
re_annotate_new_orfs <- function(orf_table) {
  orf_table <- orf_table %>%
    # Step 1: Group gene biotypes into simplified categories
    dplyr::mutate(
      gene_biotype = dplyr::case_when(
        # Unannotated transcripts from StringTie
        gene_biotype %in% c("stringtie") ~ "unannotated_gene",
        # Consolidate all pseudogene variants
        grepl("pseudogene", gene_biotype, ignore.case = TRUE) ~ "pseudogene",
        # Treat immune-related genes as protein-coding
        gene_biotype %in% c(
          "TR_V_gene", "IG_V_gene", "IG_C_gene", "TR_C_gene",
          "IG_D_gene", "IG_J_gene", "TR_D_gene", "TR_J_gene"
        ) ~ "protein_coding",
        # TEC (To be Experimentally Confirmed) as lncRNA
        gene_biotype == "TEC" ~ "lncRNA",
        # Group various small RNA biotypes together
        gene_biotype %in% c(
          "miRNA", "misc_RNA", "snRNA", "snoRNA",
          "scaRNA", "vault_RNA", "sRNA", "scRNA",
          "ribozyme"
        ) ~ "small_RNA",
        # Leave other biotypes unchanged
        TRUE ~ gene_biotype
      )
    ) %>%
    # Step 2: Refine ORF categories based on updated gene biotypes
    dplyr::mutate(
      orf_category_new = dplyr::case_when(
        # Novel ORFs in protein-coding genes become internal ORFs
        orf_category_new == "novel" &
          gene_biotype == "protein_coding" ~ "intORF",
        # Novel ORFs in pseudogenes flagged as pseudogene ORFs
        orf_category_new == "novel" &
          grepl("pseudogene", gene_biotype) ~ "pseudogene_ORF",
        # lncRNA ORFs overlapping CDS regions
        gene_biotype == "lncRNA" &
          !is.na(cds_range_similarity) ~ "ovCDS_lncRNA_ORF",
        # Other lncRNA ORFs
        gene_biotype == "lncRNA" ~ "lncRNA_ORF",
        # Remaining novel ORFs in other categories
        orf_category_new == "novel" ~ "unannotated_gene_ORF",
        # Keep existing category if no conditions match
        TRUE ~ orf_category_new
      )
    )

  # Return the updated ORF table
  return(orf_table)
}


#' Obtain p_sites from input bed file
#'
#' Reads a BED file containing P-site coordinates and filters for p0 sites,
#' returning a GRanges object of the filtered locations.
#'
#' @param column_id String. Name to assign to the ID column in the output.
#' @param bed_file String. Path to the input BED file with P0 site annotations.
#'
#' @return A GRanges object with p0 site locations and associated metadata.
process_psites <- function(column_id, bed_file) {
  # Read BED file and select relevant columns
  p0_sites <- data.table::fread(
    bed_file,
    sep = "\t",
    header = FALSE,
    select = 1:6,
    col.names = c("seqnames", "start", "stop", column_id, "codon_pos", "strand")
  ) %>%
    dplyr::filter(codon_pos == "p0")  # Keep only the p0 rows
  
  # Convert filtered data to GRanges object
  filtered_output <- p0_sites %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  return(filtered_output)
}


#' Re-annotate ORF categories based on inframe hits
#'
#' Updates ORF categories for intORFs and lncRNA_ORFs that overlap
#' known CDS P-sites, marking them as canonical or overlapping CDS,
#' and adds a boolean flag for inframe status.
#'
#' @param orf_table Data frame. Table of ORF annotations.
#' @param hit_orf_ids Character vector. ORF IDs identified as inframe hits.
#'
#' @return Updated ORF table with new category and inframe flag.
re_annotate_inframe_orfs <- function(orf_table, hit_orf_ids) {
  orf_table <- orf_table %>%
    dplyr::mutate(
      orf_category_new = dplyr::case_when(
        orf_category_new == "intORF" & orf_id %in% hit_orf_ids ~ "canonical_ORF",
        TRUE ~ orf_category_new
      )
    ) %>%
    # Add column with boolean values sindicating inframe status
    dplyr::mutate(
      Is_inframe = orf_id %in% hit_orf_ids
    )
  
  return(orf_table)
}


#' Identify and annotate inframe ORFs based on P-site overlaps
#'
#' Processes P-site BED files for reference transcripts and ORFcaller outputs,
#' finds CDS overlaps to determine inframe ORFs, and re-annotates the ORF table.
#'
#' @param orf_table Data frame: Original ORF annotation table to update.
#' @param ref_bed_file String: Path to reference transcript P-site BED file.
#' @param orfcaller_bed_file String: Path to ORFcaller P-site BED file.
#'
#' @return Updated ORF table with inframe annotation applied.
obtain_inframe_orfs <- function(orf_table, ref_bed_file, orfcaller_bed_file) {
  # Step 1: Obtain p0-sites for reference transcript CDS
  column_id <- "transcript_id"
  ref_p0 <- process_psites(column_id, ref_bed_file)
  
  # Step 2: Obtain p0-sites for ORFcaller intORFs
  column_id <- "orf_id"
  orfcaller_p0 <- process_psites(column_id, orfcaller_bed_file)
  names(orfcaller_p0) <- as.character(mcols(orfcaller_p0)[[column_id]])
  
  # Step 3: Find overlaps between reference and ORFcaller P-sites
  ol <- GenomicRanges::findOverlaps(query = ref_p0, subject = orfcaller_p0)
  
  # Step 4: Extract unique inframe ORF IDs
  hit_orf_ids <- unique(names(orfcaller_p0)[subjectHits(ol)])
  
  # Step 5: Re-annotate ORFs in the table
  orf_table <- re_annotate_inframe_orfs(orf_table, hit_orf_ids)
  
  return(orf_table)
}


# =============================================================================
# 04 | RUN ANALYSIS ----
# =============================================================================

# Step_1: Load correct ORF information depending on used ORFcaller
if (tolower(orfcaller) == "orfquant") {
  prep_orfs <- prepare_orfquant(orfquant_orfs_loc = orfs_loc)
  
} else if (tolower(orfcaller) == "price") {
  prep_orfs <- prepare_price(
    price_orfs_loc = orfs_loc,
    tx2gene = tx2gene,
    gtf_df = gtf_df
  )
} else if (tolower(orfcaller) == "ribotie") {
    prep_orfs <- prepare_ribotie(ribotie_orfs_loc = orfs_loc) 
}

# Step_2: Set same annotation style
restyled_orfs <- check_annot_style(
  orf_ranges = prep_orfs[[1]],
  annotated_gen = cds_gene
)

# Step_3: Find predicted ORF and annotated ORF location similarities
orf_cds_sim <- check_orf_cds_similarity(
  orf_ranges = restyled_orfs,
  orf_table = prep_orfs[[2]],
  annotated_gen = cds_gene,
  annotated_gen_unlist = cds_gene_unlist
)

# Step_4: Annotate ORF table
orf_table <- annotate_new_orfs(
  orf_ranges = restyled_orfs,
  orf_table = prep_orfs[[2]],
  cds_matches_grl = orf_cds_sim,
  orf_caller = orfcaller
)

# Step_5: Rename gene biotypes and orf category
orf_table <- re_annotate_new_orfs(orf_table = orf_table)

# Step_6: Intorf re-annotation
orf_table <- obtain_inframe_orfs(orf_table, ref_bed_file, orfcaller_bed_file)


# =============================================================================
# 05 | SAVE OUTPUT ----
#   * Write final ORF table to csv file
# =============================================================================

# Write results to table
write.table(
  orf_table,
  file = file.path(".", paste0(orfcaller, "_orfs.csv")),
  quote = FALSE,
  row.names = FALSE,
  sep = ","
)

