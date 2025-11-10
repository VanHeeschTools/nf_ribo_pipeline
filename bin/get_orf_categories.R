#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(rtracklayer)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(Biostrings)
    library(stringr)
})

# Obtain input arguments
args <- commandArgs(trailingOnly=TRUE)

gtf_file <- args[1]
orf_file <- args[2]
cds_orf_bed_file <- args[3]
cds_id_gene_file <- args[4] 
bsgenome_path <- args[5]
orfcaller <- args[6]

# Create txdb file using reference gtf
txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)

## Define functions

#' Extract lncRNA transcript IDs from a GTF file
#'
#' This function imports a GTF file, filters for transcript entries
#' annotated as lncRNAs, and returns a data frame of corresponding transcript IDs.
#'
#' @param gtf_file Path to a GTF annotation file.
#'
#' @return A data frame with transcript id and bool that equals TRUE
extract_transcript_biotypes <- function(gtf_file) {
    # Import GTF annotation
    gtf <- rtracklayer::import(gtf_file)
    
    # Keep only transcript entries
    tx <- gtf[gtf$type == "transcript"]
    
    # Create a data frame with transcript ID and transcript biotype
    tx_df <- data.frame(
        tx_id = tx$transcript_id,
        biotype = if (!is.null(tx$transcript_biotype)) {
        tx$transcript_biotype
        } else if (!is.null(tx$gene_biotype)) {
        tx$gene_biotype # Set gene_biotype if transcript_biotype not found
        } else {
        NA_character_
        },
        stringsAsFactors = FALSE
    )
    
    tx_df <- tx_df %>% dplyr::mutate(
        biotype = dplyr::case_when(
        # Unannotated transcripts from StringTie
        biotype %in% c("stringtie") ~ "unannotated_gene",
        # Consolidate all pseudogene variants
        grepl("pseudogene", biotype, ignore.case = TRUE) ~ "pseudogene",
        # Treat immune-related genes as protein-coding
        biotype %in% c(
            "TR_V_gene", "IG_V_gene", "IG_C_gene", "TR_C_gene",
            "IG_D_gene", "IG_J_gene", "TR_D_gene", "TR_J_gene"
        ) ~ "protein_coding",
        # TEC (To be Experimentally Confirmed) as lncRNA
        biotype == "TEC" ~ "lncRNA",
        # Group various small RNA biotypes together
        biotype %in% c(
            "miRNA", "misc_RNA", "snRNA", "snoRNA",
            "scaRNA", "vault_RNA", "sRNA", "scRNA",
            "ribozyme"
        ) ~ "small_RNA",
        # Leave other biotypes unchanged
        TRUE ~ biotype
        )
    )
    
    return(tx_df)
}

#' Obtain in-frame ORFs overlapping CDS regions using P0 sites
#'
#' This function reads the P0 BED file between CDS and ORF regions and returns 
#' the unique pairs of `cds_id` and `orf_id` that have an in-frame overlap.
#'
#' @param cds_orf_bed_file Path to the BED file containing CDS–ORF P0 
#'  intersections.
#'
#' @return A dataframe with two columns: cds_id and orf_id, containing the 
#'  distinct in-frame overlaps.
obtain_inframe_orfs <- function(cds_orf_bed_file) {
    # Read BED intersection file and set column names
    cds_orf_intersect <- data.table::fread(
        cds_orf_bed_file,
        col.names = c(
        "chr", "start", "end", "cds_id", "frame", "strand", "nt_position",
        "chr_match", "start_match", "end_match", "orf_id", "frame_orf",
        "strand_orf", "nt_position_orf"
        )
    ) %>%
    # Only keep locations where reference and ORF have a p0 location
    filter(frame == "p0" & frame_orf == "p0") %>% 
    # Keep unique CDS ORF p0 pairs (remove duplicates)
    dplyr::distinct(cds_id, orf_id)
    
    return(cds_orf_intersect)
}

#' Load ORFcaller gtf file
#'
#' This function loads the ORFcaller gtf file and manipulates the data to obtain
#' required columns
#'
#' @param gtf_file_path Path to the ORFcaller GTF file.
#'
#' @return A gtf dataframe 
#'
load_orfcaller_gtf <- function(gtf_file, orf_gtf_file, txdb, orfcaller){
    
    orf_gtf_df <- rtracklayer::import(orf_gtf_file) %>%
        as.data.frame()
    
    # Summarize CDS entries per ORF
    orf_list <- orf_gtf_df %>%
        dplyr::filter(type == "CDS") %>%             # Keep only CDS features
        dplyr::group_by(ORF_id) %>%                  # Group by ORF identifier
        dplyr::summarise(
        chrm   = unique(as.character(seqnames))[1],  # Extract chromosome
        strand = unique(as.character(strand))[1],    # Extract strand
        starts = paste(sort(start), collapse = ","), # Combine sorted starts
        ends   = paste(sort(end), collapse = ","),    # Combine sorted ends
        # Get ORFcaller given transcript_id to later obtain gene metadata
        transcript_id = unique(as.character(transcript_id))[1]
        ) %>%
        dplyr::ungroup() %>%
        mutate(
        orf_id = ORF_id,
        orf_start = str_extract(starts, "^[^,]+") %>% as.numeric(),
        orf_end  = str_extract(ends, "[^,]+$") %>% as.numeric(),
        summary_id = paste0(orfcaller, "_", chrm, ":", orf_start, "-", orf_end, "_", strand), 
        orfcaller = orfcaller
        )%>%
        rowwise() %>%
        mutate(
        orf_width = sum(as.numeric(str_split(ends, ",")[[1]]) -
                            as.numeric(str_split(starts, ",")[[1]]) + 1)
        ) %>%
        dplyr::select(-ORF_id) %>%
        as.data.frame()

    # Load gene_meta_data to add gene_id, gene_name and gene_biotype
    reference_gtf_df <- rtracklayer::import(gtf_file) %>%
    as.data.frame() %>%
    mutate(
        transcript_id = case_when(
        orfcaller == "RiboTIE" & str_starts(transcript_id, "TCONS_") ~ 
            str_extract(transcript_id, "^[^_]+_[^_]+"),      # In case of a TCONS transcript keep everything untill second underscore
        orfcaller == "RiboTIE" ~ 
            str_replace(transcript_id, "_.*$", ""),          # Otherwise only keep untill the first underscore
        TRUE ~ ranscript_id                                  # If not a RiboTIE ORF leave unchanged
        )
    )

    gene_meta <- reference_gtf_df %>%
        as.data.frame() %>%
        dplyr::filter(type == "transcript") %>%
        dplyr::select(transcript_id, gene_id, gene_name, gene_biotype) %>%
        dplyr::distinct() 
    
    orf_list <- orf_list %>%
        dplyr::left_join(gene_meta, by = "transcript_id") %>% 
        dplyr::select(-transcript_id)
    
    
    return(orf_list)
}

#' Match ORFs to annotated transcripts and measure overlap
#'
#' This function matches ORFs from an ORFcaller output to transcript models
#' from a GTF file and calculates the measured overlap width for each
#' ORF–transcript pair. All transcripts with which the ORF has overlap will be 
#' found.
#'
#' @param orf_file Path to the ORFcaller-derived GTF file.
#' @param gtf_file Path to the reference GTF annotation file.
#'
#' @return A dataframe with one row per ORF–transcript match
#'
match_orfs_to_transcripts <- function(orf_file, gtf_file, txdb) {
    # Get all transcripts from GTF
    transcript_list <- exonsBy(txdb, by = "tx", use.names = TRUE)
    
    orf_list_expanded <- orf_list %>%
        dplyr::select(orf_id, chrm, strand, starts, ends) %>%
        separate_rows(starts, ends) %>%
        mutate(starts = as.numeric(starts),
            ends = as.numeric(ends))
    
    # Convert ORFs to GRanges
    orf_grl <- GRanges(
        seqnames = orf_list_expanded$chrm,
        ranges   = IRanges(start = orf_list_expanded$starts,
                        end = orf_list_expanded$ends),
        strand   = orf_list_expanded$strand) %>%
        setNames(orf_list_expanded$orf_id) %>%
        split(orf_list_expanded$orf_id)
    
    # Find ORFs that overlap completely with transcript
    hits <- findOverlaps(orf_grl, transcript_list, type = "within") %>%
        as.data.frame() %>%
        right_join(data.frame(orf_id = names(orf_grl)) %>%
                    mutate(queryHits = row_number()), "queryHits") %>%
        left_join(data.frame(tx_id = names(transcript_list)) %>%
                    mutate(subjectHits = row_number()), "subjectHits") %>%
        left_join(orf_list %>% dplyr::select(orf_id, chrm, orf_start, orf_end, orf_width,  
                                            strand), "orf_id") %>%
        mutate(tx_orf_id = paste0(tx_id, "__", orf_id)) %>%
        # Filter out ncORFs located on novel transcripts
        dplyr::filter(str_starts(tx_id, "ENST"))
    
    # Determine width for each ORF-transcript match (to check for correct
    # transcript structure)
    hits_with_range <- GRanges(
        seqnames = hits$chrm,
        ranges   = IRanges(start = hits$orf_start,
                        end = hits$orf_end),
        orf_id = hits$orf_id)%>%
        pintersect(transcript_list[hits$tx_id], drop.nohit.ranges = TRUE) %>%
        setNames(hits$tx_orf_id) %>%
        data.frame() %>%
        group_by(tx_orf_id = group_name) %>%
        summarise(measured_width = sum(width)) %>%
        inner_join(hits, "tx_orf_id")
    
    return(hits_with_range)
}

#' Summarize CDS coordinate limits per transcript
#'
#' This function reads an RDS file containing CDS annotations, groups them
#' by transcript ID, and determines the start and end for each transcript.
#'
#' @param cds_id_gene_file Path to an RDS file containing CDS information.
#'
#' @return A tibble with one row per transcript containing:
#'     Transcript ID, minimum CDS start coordinate, maximum CDS end coordinate
#'
summarize_cds_tx_limits <- function(cds_id_gene_file) {
    # Load CDS annotation table from RDS file
    cds_data <- readRDS(cds_id_gene_file)
    
    # Summarize CDS coordinate limits (min start, max end) per transcript
    cds_id_tx_limits <- cds_data %>%
        dplyr::group_by(tx_id) %>%
        dplyr::summarise(
        cds_start = min(start),
        cds_end   = max(end),
        .groups = "drop"
        )
    
    return(cds_id_tx_limits)
}

#' Categorize ORFs based on CDS overlap type
#'
#' This function integrates ORF–CDS intersection data with transcript
#' overlap metrics and CDS information, assigning each ORF to a category:
#'
#' @param cds_orf_intersect Data frame containing unique CDS–ORF intersections
#' @param hits_with_range Data frame with ORF–transcript matches and measured
#'   overlap widths
#' @param cds_id_tx_limits Data frame with CDS start and end coordinates per
#'   transcript
#'
#' @return A tibble with one row per ORF, containing ORF_id and category
#'
categorize_orf_cds_overlap <- function(cds_orf_intersect,
                                        hits_with_range,
                                        cds_id_tx_limits) {

    # Determine category for all transcripts for which the ORF has exon overlap or p0 CDS overlap
    cds_overlap_orfs <- cds_orf_intersect %>%
        # Align column names for joining
        dplyr::rename(tx_id = cds_id) %>%
        # Join overlap and coordinate summary data
        dplyr::left_join(hits_with_range, by = c("orf_id", "tx_id")) %>%
        dplyr::left_join(cds_id_tx_limits, by = "tx_id") %>%
        # Classify overlap category
        dplyr::mutate(
        cds_overlap_cat = dplyr::case_when(
            # If there are no matching transcripts found its an incomplete overlap
            is.na(tx_orf_id) ~ "Incomplete_overlap",
            # It is annotated if size and coordinates match completely
            measured_width == orf_width & cds_start == orf_start & cds_end == orf_end ~ 
            "ORF-annotated",
            # If size is equal but not all coordinates its an NC-variant
            measured_width == orf_width ~ "NC-variant",
            # If width is not equal it means part of the CDS is missing so it is an CDS-isoform
            measured_width != orf_width ~ "CDS-isoform") ) %>%  
        # Keep only one category per ORF
        dplyr::select(orf_id, tx_id, cds_overlap_cat)

    # Obtain ORFs that have overlap with kown p0 sites  
    cds_overlap_orfs_summarized <- cds_overlap_orfs %>% 
        mutate(
        cds_overlap_cat = factor(
            cds_overlap_cat,
            levels = c("ORF-annotated", "NC-variant", "CDS-isoform", "Incomplete_overlap")
        )
        ) %>% 
        arrange(orf_id, cds_overlap_cat) %>%
        group_by(orf_id) %>%
        summarise(has_cds_overlap = dplyr::first(cds_overlap_cat), .groups = "drop")
    
    return(list(
        cds_overlap_orfs = cds_overlap_orfs,
        cds_overlap_orfs_summarized = cds_overlap_orfs_summarized
    ))
}

#' Classify ORFs by transcript overlap
#'
#' This function classifies ORFs based on how they overlap with known transcripts,
#' CDS, and transcript biotypes.
#'
#' @param orf_list Data frame of ORFs with coordinates and strand information.
#' @param hits_with_range Data frame of ORF–transcript overlaps, including measured widths.
#' @param cds_overlap_orfs Data frame mapping ORFs to CDS overlap categories.
#' @param cds_id_tx_limits Data frame with per-transcript CDS start and end coordinates.
#' @param lncRNA_ids Data frame with transcript IDs annotated as lncRNAs.
#'
#' @return A dataframe containing the ORFs and their given categories
#'
classify_orfs <- function(orf_list, hits_with_range, cds_overlap_orfs, cds_overlap_orfs_summarized, cds_id_tx_limits) {
    
    
    # Obtain all ORFs with transcript on which it could be an isoform
    orf_isoform_transcripts <- cds_overlap_orfs %>% 
        dplyr::filter(cds_overlap_cat=="CDS-isoform") %>% 
        dplyr::group_by(orf_id) %>% 
        dplyr::summarise(tx_id = paste0(tx_id, collapse = "__")) %>% 
        dplyr::rename(cds_isoform_transcripts = tx_id)
    
    
    orf_list_cat <- orf_list %>%
        ## Join in transcript matches with matching measured width
        dplyr::left_join(
        hits_with_range %>%
            dplyr::filter(measured_width == orf_width) %>%
            dplyr::select(orf_id, tx_id),
        by = "orf_id"
        ) %>%
        ## Add CDS overlap category and transcript boundaries
        dplyr::left_join(cds_overlap_orfs, by = c("tx_id", "orf_id")) %>%
        dplyr::left_join(cds_id_tx_limits, by = "tx_id") %>%
        dplyr::left_join(cds_overlap_orfs_summarized, by = "orf_id") %>% 
        dplyr::left_join(transcript_biotype_ids, by = "tx_id") %>%

        # Change cds_overlap column to show if it has p0 overlap or not
        dplyr::rename(has_p0_cds_overlap = has_cds_overlap) %>%
        dplyr::mutate(has_p0_cds_overlap = !is.na(has_p0_cds_overlap))%>%

        ## Assign ORF category
        dplyr::mutate(
        
        orf_cat = dplyr::case_when(
            # No transcript assignment
            is.na(tx_id) ~ "undefined",
            
            # Has both cds_overlap_cat and lncRNA biotype
            has_p0_cds_overlap == "Has_p0_CDS_overlap" & biotype == "lncRNA" ~ "ovCDS_lncRNA-ORF",
            #TODO: what about out of frame CDS overlap
            
            #Check if novel-ORFs can be intORFs
            biotype == "lncRNA" ~ "lncRNA-ORF",
            biotype == "pseudogene" ~ "pseudogene-ORF",
            
            # Already has CDS overlap category
            !is.na(cds_overlap_cat) ~ cds_overlap_cat,
            
            # Novel transcript ORFs (e.g. from StringTie)
            stringr::str_starts(tx_id, "TCONS") ~ "novel-ORF",
            
            # Processed transcript without CDS
            is.na(cds_start) ~ "Processed_transcript_ORF",
            
            strand == "+" & orf_end < cds_start ~ "uORF",
            strand == "+" & orf_start < cds_start & orf_end >= cds_start ~ "uoORF",
            strand == "+" & orf_start > cds_start & orf_end < cds_end ~ "intORF", # Will only be found if not inframe
            strand == "+" & orf_start <= cds_end & orf_end > cds_end ~ "doORF",
            strand == "+" & orf_start > cds_end ~ "dORF",
            
            ## Minus-strand ORFs relative to CDS 
            strand == "-" & orf_end < cds_start ~ "dORF",
            strand == "-" & orf_start < cds_start & orf_end >= cds_start ~ "doORF",
            strand == "-" & orf_start > cds_start & orf_end < cds_end ~ "intORF", # Will only be found if not inframe
            strand == "-" & orf_start <= cds_end & orf_end > cds_end ~ "uoORF",
            strand == "-" & orf_start > cds_end ~ "uORF",
            
            # Fallback if no other condition matched
            .default = "unassigned"),
        # Add factor for determining category in case of multiple transcripts
        orf_cat = factor(orf_cat, c("ORF-annotated", "NC-variant", "ovCDS_lncRNA-ORF",
                                    "CDS-isoform", "lncRNA-ORF","pseudogene-ORF",
                                    "Incomplete_overlap", "intORF", "uoORF",
                                    "doORF", "uORF", "dORF",
                                    "Processed_transcript_ORF", "novel-ORF",
                                    "undefined", "unassigned")))
    
    # Summarise orf_table to get one row per ORF 
    orf_list_sum <- orf_list_cat %>%
        arrange(orf_id, orf_cat) %>%
        group_by(orf_id, summary_id, chrm,orf_start, orf_end, strand, starts, ends, 
        has_p0_cds_overlap,orfcaller, gene_id, gene_name, gene_biotype ) %>%
        summarise(tx_id = paste0(tx_id, collapse = "__"),
                transcript_biotype_all = paste0(biotype, collapse = "__"),
                orf_biotypes_all = paste0(orf_cat, collapse = "__"),
                orf_biotype_single = dplyr::first(orf_cat),
                .groups = "drop")%>%
        mutate(
        summary_id = paste0(summary_id, "_", orf_biotype_single,"_",orf_id),
        )%>% 
        dplyr::left_join(orf_isoform_transcripts, by = "orf_id") %>% 
        dplyr::mutate(cds_isoform_transcripts = if_else(
            is.na(cds_isoform_transcripts), 
            "None", cds_isoform_transcripts))
}

load_bsgenome_library <- function(bsgenome_path){
    # Add the directory to the library search path
    paths <- c(bsgenome_path, .libPaths())
    .libPaths(paths)
    
    # Automatically detect the BSgenome package name in that directory
    bsgenome_dirs <- list.dirs(bsgenome_path, full.names = FALSE, 
                                recursive = FALSE)
    bsgenome_name <- bsgenome_dirs[grepl("^BSgenome\\.", bsgenome_dirs)]
    # Function to extract the BSgenome object from the package namespace
    load_bsgenome <- function(pkg) {
        ns <- asNamespace(pkg)
        bsgenome_obj <- Filter(function(x) inherits(get(x, envir = ns), "BSgenome"),
                            ls(ns))
        get(bsgenome_obj[[1]], envir = ns)
    }
    
    # Load the genome object
    genome <- load_bsgenome(bsgenome_name)
    return(genome)
    
}

#' Generate ORF sequences and protein translations
#'
#' This function takes an ORF table and a reference genome, merges exon sequences
#' per ORF strand-aware, translates them to protein sequences, and returns a data frame
#' with nucleotide and protein sequences including start and stop codons.
#'
#' @param orf_table A data frame containing the ORF information.
#' @param genome A `BSgenome` object representing the reference genome.
#' 
#' @return The orf_table dataframe plus added Proteins sequence and start, stop
#' codon info
#'   
generate_orf_sequences <- function(orf_table, genome) {
    
    # Expand ORF exons and order them strand-aware
    orf_long <- orf_table %>%
        separate_rows(starts, ends, sep = ",") %>%
        mutate(
        start = as.integer(starts),
        end   = as.integer(ends)
        ) %>%
        arrange(orf_id, if_else(strand == "+", start, -start))
    
    # Build GRanges object from ORF exon coordinates
    orf_gr <- makeGRangesFromDataFrame(
        orf_long,
        seqnames.field = "chrm",
        start.field    = "start",
        end.field      = "end",
        strand.field   = "strand",
        keep.extra.columns = TRUE
    )
    
    # Split GRanges by ORF
    orf_grl <- split(orf_gr, orf_gr$orf_id)
    
    # Extract nucleotide sequences for all exons
    dna_list <- getSeq(genome, orf_grl)  # DNAStringSetList
    
    # Collapse exon sequences per ORF into single DNAStringSet
    df_seq <- data.frame(
        orf_id = rep(names(dna_list), lengths(dna_list)),
        seq    = as.character(unlist(dna_list)),
        stringsAsFactors = FALSE
    )
    
    dna_seqs <- DNAStringSet(
        tapply(df_seq$seq, df_seq$orf_id, paste0, collapse = "")
    )
    names(dna_seqs) <- names(tapply(df_seq$seq, df_seq$orf_id, paste0, collapse = ""))
    
    # Translate to protein sequences
    protein_aa <- Biostrings::translate(dna_seqs)
    
    # Convert to plain character
    protein_seqs <- as.character(protein_aa)
    
    # Remove trailing stop codons (*) if any
    protein_seqs <- sub("\\*$", "", protein_seqs)  
    # Build final ORF table
    orf_final <- tibble(
        orf_id      = names(dna_seqs),
        dna_seq     = as.character(dna_seqs),
        protein_seq = protein_seqs,
        protein_length = length(protein_seq),
        start_codon = substr(dna_seq, 1, 3),
        stop_codon  = substr(dna_seq, nchar(dna_seq)-2, nchar(dna_seq))
    ) %>%
        left_join(distinct(orf_table, orf_id, .keep_all = TRUE), by = "orf_id") %>%
        dplyr::mutate(protein_length = nchar(protein_seq)) %>%
        dplyr::select(orf_id, summary_id, gene_id, gene_name, gene_biotype, 
                    protein_seq, protein_length, dna_seq, start_codon, stop_codon, 
                    chrm, orf_start, orf_end, strand, starts, ends, 
                    has_p0_cds_overlap, cds_isoform_transcripts, orfcaller,
                    tx_id, transcript_biotype_all, orf_biotypes_all, 
                    orf_biotype_single) %>%
        dplyr::mutate(starts = gsub(",", "_", starts)) %>%
        dplyr::mutate(ends = gsub(",", "_", ends)) %>%
        dplyr::arrange(chrm, orf_start, orf_end) # Sort on genomic coordinates
    
    return(orf_final)
    
}

# Run functions:

# Extract the biotypes of the reference transcripts
transcript_biotype_ids <- extract_transcript_biotypes(gtf_file)

# Load p0 overlap bed file to determine ORF frame
cds_orf_intersect <- obtain_inframe_orfs(cds_orf_bed_file)

# Load ORFcaller ORFS
orf_list <- load_orfcaller_gtf(gtf_file, orf_file, txdb, orfcaller)

# Obtains overlap of ORF CDS with reference exons
hits_with_range <- match_orfs_to_transcripts(orf_file, gtf_file, txdb)

# Obtains starts and stops of reference CDS from RDS created in p0 bed file script
cds_id_tx_limits <- summarize_cds_tx_limits(cds_id_gene_file)

# Classify ORF CDS on overlap with reference exons
cds_overlap_list <- categorize_orf_cds_overlap(
    cds_orf_intersect,
    hits_with_range,
    cds_id_tx_limits)

cds_overlap_orfs <- cds_overlap_list$cds_overlap_orfs
cds_overlap_orfs_summarized <- cds_overlap_list$cds_overlap_orfs_summarized

# Classify remaning ORFs using transcript biotype and reference CDS
orf_table <- classify_orfs(orf_list, hits_with_range, 
                            cds_overlap_orfs, cds_overlap_orfs_summarized,
                            cds_id_tx_limits)

# Load BSgenome library to obtain ORF nucleotide sequence
genome <- load_bsgenome_library(bsgenome_path)

# Add Protein sequence to ORF table and finalise table
orf_results <- generate_orf_sequences(orf_table, genome)

# Write results to table
write.table(
    orf_results,
    file = file.path(".", paste0(orfcaller, "_orfs.csv")),
    quote = FALSE,
    row.names = FALSE,
    sep = ","
)