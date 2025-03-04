#!/usr/bin/env Rscript

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
})

# This script creates new annotations for PRICE and ORFquant ORFs based on the annotated
# reference ORFs

args <- commandArgs(trailingOnly = TRUE)
orfs_loc <- args[1]
gtf <- args[2]
orfcaller <- args[3]
annotation_provider <- args[4]
gencode_uniprot_file <- args[5]
uniprot_protein_fasta_loc <- args[6]
package_install_loc <- args[7]
bsgenome_path <- args[8]
cpus <- args[9]


# Load custom BSgenome package
paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

load_bsgenome <- function(pkg) {
  ns <- asNamespace(pkg)
  bsgenome_obj <- Filter(function(x) inherits(get(x, envir = ns), "BSgenome"), ls(ns))
  get(bsgenome_obj[[1]], envir = ns)
}
genome <- load_bsgenome(basename(bsgenome_path))

# Extract CDS regions of annotated genes from txdb
txdb <- makeTxDbFromGFF(gtf, format = "gtf")
#txdb <- AnnotationDbi::loadDb(txdb_loc)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
cds_gene <- GenomicFeatures::cdsBy(txdb, "gene")
cds_gene_unlist <- unlist(cds_gene)


# FUNCTIONS
prepare_price <- function(price_orfs_loc, tx2gene) {
  #price_orfs_loc<-orfs_loc
  # This functions imports the PRICE ORFs into R
  # Returns a list with the ORF genomci ranges (orf_ranges)
  # and the metadata info per ORF (price_orf_df)
  
  # Import PRICE ORF definitions
  price_orfs <- rtracklayer::import.bed(price_orfs_loc)
  orf_ranges <- rtracklayer::blocks(price_orfs) # uses BED file blocks to establish exons
  
  # Fix gene IDs for transcript IDs (which are correct)
  price_orf_df <- as.data.frame(price_orfs) %>%
    dplyr::select(seqnames,start,end,width,strand,name) %>%
    dplyr::mutate(gene_id = stringr::str_split_i(name, "__", i = 1),
                  transcript_id = stringr::str_split_i(name, "__", i = 2),
                  start_codon = stringr::str_split_i(name, "__", i = 5)) %>%
    dplyr::left_join(tx2gene, by = c("transcript_id" = "TXNAME")) %>%
    dplyr::select(!(gene_id))%>%
    dplyr::rename(gene_id = GENEID)

  
  # Obtain orf AA sequences for Price orfs
  # Check if genetic code is needed here and what it does
  # Create a custom genetic code (if required)
  GENETIC_CODE_PRICE <- Biostrings::GENETIC_CODE
  attr(GENETIC_CODE_PRICE, "alt_init_codons") <- unique(price_orf_df$start_codon)
  
  # Split the ORF table by strand for simplicity
  price_plus <- price_orf_df$`name`[price_orf_df$strand == "+"]
  price_minus <- price_orf_df$`name`[price_orf_df$strand == "-"]
  
  # Sort the minus strand (if necessary)
  orf_minus_sorted <- GenomicRanges::sort(orf_ranges[which(names(orf_ranges) %in% price_minus)], decreasing = TRUE)
  
  # Define a function to translate sequences and return protein sequences
  get_protein_seq <- function(orf_names, genome) {
    lapply(Biostrings::getSeq(x = genome, names = orf_ranges[which(names(orf_ranges) %in% orf_names),]), 
           function(seq) {
             Biostrings::toString(Biostrings::translate(Biostrings::DNAString(paste(seq, collapse = "")), genetic.code = GENETIC_CODE_PRICE))
           })
  }
  
  # Get protein sequences for both strands
  price_protein_plus <- get_protein_seq(price_plus, genome)
  price_protein_minus <- get_protein_seq(names(orf_minus_sorted), genome)
  
  # Combine results into a data frame
  price_proteins <- data.frame(
    name = c(names(price_protein_plus), names(price_protein_minus)),
    Protein = c(unlist(price_protein_plus), unlist(price_protein_minus))
  )
  
  # Join with the original dataframe
  price_orf_df <- price_orf_df %>% left_join(price_proteins, by = "name")
  
  return(list(orf_ranges,
              price_orf_df))
}

prepare_orfquant <- function(orfquant_orfs_loc) {
  
  # Imports the ORFquant object and outputs a list of ORF ranges (orf_ranges), 
  # and the ORF metadata (orfs_table) linked to the ORF ranges
  
  # Load ORFquant file and select ORF definitions
  orfquant_orfs <- get(load(orfquant_orfs_loc))
  orfs_tx_df <-  data.frame(S4Vectors::mcols(orfquant_orfs$ORFs_tx)[, c("ORF_id_tr", "Protein", "gene_id", "gene_biotype", "gene_name", "transcript_id", "transcript_biotype", "ORF_category_Tx", "ORF_category_Gen", "P_sites_raw", "P_sites_raw_uniq")]) %>%
    dplyr::distinct()
  
  # Generate table of genomic ORF locations
  orf_table <- data.frame(orfquant_orfs$ORFs_gen) %>%
    dplyr::mutate(ORF_id_tr = names(orfquant_orfs$ORFs_gen)) %>%
    dplyr::group_by(ORF_id_tr) %>%
    dplyr::mutate(ORF_ranges = paste0(seqnames, ":", min(start), "-", max(end))) %>%
    dplyr::select(c("ORF_id_tr", "ORF_ranges")) %>%
    dplyr::distinct() %>%
    dplyr::left_join(orfs_tx_df, by = c("ORF_id_tr")) # Add column for uniprot
  
  orf_ranges <- split(orfquant_orfs$ORFs_gen, names(orfquant_orfs$ORFs_gen))
  
  return(list(orf_ranges,
              orf_table
  ))
}

check_annot_style <- function(orf_ranges, annotated_gen) {
  
  # Convert the seqlevels style of a GRange so that they can be compared with one
  # another. Only checks for NCBI and Ensembl.
  # Input: orf_ranges and CDS definitions of the ORF, only checks the seqlevels of
  #        these GRanges
  # Output: orf_ranges with updated seqlevels
  
  # This is hard-coded, I would recommend checking the seqlevels of the
  # Called ORFs and use those as a basis for the conversion
  ensembl_seqlevels <- c("1","2","3","4","5","6","7","8","9","10","11","12",
                         "13","14","15","16","17","18","19","20","21","22","X")
  ncbi_seqlevels <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                      "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                      "chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                      "chrX")
  # Required to reannotate the seqlevels
  names(ensembl_seqlevels) <- ncbi_seqlevels
  names(ncbi_seqlevels) <- ensembl_seqlevels
  
  # Checks whether to change the PRICE seqnames to NCBI style instead of Ensembl
  check_cds <- ifelse(any(grepl("chr",
                                GenomeInfoDb::seqlevels(annotated_gen[1]))),
                      "ncbi",
                      "ensembl")
  check_caller <- ifelse(any(grepl("chr",
                                   GenomeInfoDb::seqlevels(orf_ranges[1]))),
                         "ncbi",
                         "ensembl")
  
  if(!(check_cds == check_caller)) {
    print(paste("switch to",check_cds))
    if(check_cds == "ncbi") {
      GenomeInfoDb::seqlevels(orf_ranges) <- ncbi_seqlevels
    } else if (check_cds == "ensembl") {
      GenomeInfoDb::seqlevels(orf_ranges) <- ensembl_seqlevels
    }
  }
  return(orf_ranges)
}

check_orf_cds_similarity <- function(orf_ranges, orf_table, annotated_gen, annotated_gen_unlist) {
  
  # This function 
  # Input: orf_ranges (Granges of orfcaller ORFs)
  #        orf_table (orf metadata)
  #        annotated_gen (Granges of annotated ORFs)
  #        annotated_gen_unlist (Grangeslist of annotated ORFs)
  # Output: result_list (list of )
  
  # Find overlaps between called ORF and annotated ORF
  overlaps <- GenomicRanges::findOverlaps(orf_ranges, annotated_gen)
  
  # Calculate overlap between the two ORF types
  overlap_width <- sum(width(GenomicRanges::intersect(orf_ranges[S4Vectors::queryHits(overlaps)], 
                                                      annotated_gen[S4Vectors::subjectHits(overlaps)])))
  
  # Save in DF
  overlap_df <- data.frame(queryIdx = S4Vectors::queryHits(overlaps), 
                           subjectIdx = S4Vectors::subjectHits(overlaps),
                           overlapWidth = overlap_width)
  
  # Only keep the ORFs with the highest overlap
  max_overlaps <- overlap_df[order(overlap_df$queryIdx, -overlap_df$overlapWidth),]
  max_overlaps <- max_overlaps[!duplicated(max_overlaps$queryIdx),]
  
  query_idx <- max_overlaps$queryIdx
  subject_idx <- max_overlaps$subjectIdx
  
  # Populate new DF with filtered overlaps
  selected_overlaps <- data.frame(
    queryHits = 1:length(orf_ranges),
    subjectHits = rep(NA, length(orf_ranges))
  )
  
  selected_overlaps$subjectHits[selected_overlaps$queryHits %in% query_idx] <- subject_idx
  
  result_list <- GenomicRanges::GRangesList(rep(list(GenomicRanges::GRanges()), length(orf_ranges)))
  names(result_list) <- names(orf_ranges)
  
  # Annotate the ORFs with overlap correctly
  non_na_indices <- !is.na(selected_overlaps$subjectHits)
  result_list[selected_overlaps$queryHits[non_na_indices]] <- annotated_gen[selected_overlaps$subjectHits[non_na_indices]]
  no_overlap_idx <- lengths(result_list) == 0
  no_overlap_names <- names(which(no_overlap_idx))
  
  # Annotate the ORFs with no CDS overlap
  result_list[no_overlap_idx] <- GenomicRanges::GRangesList(lapply(no_overlap_names, function(name) {
    # You need orf_table here, which contains mappings between ORF IDs and parent gene IDs
    orf_parent_gene <- orf_table$gene_id[match(name, orf_table$ORF_id_tr)] 
    # Turns out you dont need to find nearest CDS regions using `nearest()`, could just use the 
    # parent gene ID -> an ORF can't be a dORF or uORF if it's in a different gene
    cds_parent_gene <- annotated_gen_unlist[which(names(annotated_gen_unlist) == orf_parent_gene)]
    return(cds_parent_gene)
  }))
  
  return(result_list)
  
}

annotate_new_orfs <- function(orf_ranges, orf_table, cds_matches_grl, orf_caller) {
  # Compare CDS vs new ORF, and make a new annotation based on a virtual longest
  # possible ORF of the annotated set.  
  # Input: orf_ranges (Granges of orf_caller ORFs)
  #        orf_table (metadata of orf_caller ORFs)
  #        cds_matches_grl ()
  #        orf_caller (This is important as ORFquant does not include the STOP codon in the sequence)
  # Output: orf_table with new columns showing start codon, orf coord similarity and the new
  #         orf annotation type.
  
  # ORFquant does not include stops in the ORF, so we trim the CDS matches
  
 
  # Convert ORF caller name to lowercase for consistency
  orf_caller <- tolower(orf_caller)
  adjust_stop <- ifelse(orf_caller == "orfquant", 3, 0)  # Adjust stop codon handling based on ORF caller
  
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
  classify_orf <- function(start, stop, ann_start, ann_stop) {
    case_when(
      stop == ann_stop & start == ann_start ~ "ORF_annotated",
      stop == ann_stop & start < ann_start ~ "N_extension",
      stop == ann_stop & start > ann_start ~ "N_truncation",
      stop != ann_stop & start < ann_start & stop < ann_stop ~ "overl_uORF",
      stop != ann_stop & start < ann_start & stop < ann_start ~ "uORF",
      stop != ann_stop & start < ann_start & stop > ann_stop ~ "NC_extension",
      stop != ann_stop & start > ann_start & stop > ann_stop ~ "overl_dORF",
      stop != ann_stop & start > ann_stop & stop > ann_stop ~ "dORF",
      stop != ann_stop & start > ann_start & stop < ann_stop ~ "nested_ORF",
      stop != ann_stop & start == ann_start & stop < ann_stop ~ "C_truncation",
      stop != ann_stop & start == ann_start & stop > ann_stop ~ "C_extension"
    )
  }
  
  # Apply classification function
  orf_category <- classify_orf(orf_start, orf_stop, ann_start, ann_stop)
  
  # Mark novel ORFs if they do not match any CDS
  orf_category[lengths(cds_matches_grl) == 0] <- "novel"
  cds_range_similarity[lengths(cds_matches_grl) == 0] <- NA
  
  # Extract start codon, specific handling for price
  start_codon <- ifelse(orf_caller == "price", stringr::str_split_i(names(orf_ranges), "__", i = 5), "ATG")
  
  # Create new ORF category data frame
  new_category_df <- data.frame(
    orf_id = names(orf_ranges),
    start_dif = abs(orf_start - ann_start),
    cds_range_similarity = as.numeric(cds_range_similarity),
    stop_same = orf_stop == ann_stop,
    orf_category_new = orf_category,
    start_codon = ifelse(tolower(orf_caller) == "price",
                         stringr::str_split_i(names(orf_ranges), "__", i = 5),
                         "ATG")
  ) %>%
    dplyr::mutate(
      orf_cds = cds_range_similarity >= 0.9 & cds_range_similarity <= 1.1,
      start_check = start_dif > 0,
      same_as_cds = start_dif < 99 & start_check & stop_same & orf_cds &
        (tolower(orf_caller) != "price" | start_codon != "ATG")
    )
  
  # Remove start_codon column for "price" ORF caller
  if (orf_caller == "orfquant"){
    orf_table <- orf_table %>%
      dplyr::left_join(new_category_df, by = c("ORF_id_tr" = "orf_id"))
  }
  
  if (orf_caller == "price") {
    orf_table <- orf_table %>%
      dplyr::select(!start_codon) %>%
      dplyr::left_join(new_category_df, by = c("name" = "orf_id"))
  }
  
  return(orf_table)
}

annotate_uniprot_id <- function(annotation_provider, gencode_uniprot_file, orf_table) {  
  # Annotate ORFs with UniProt ID where possible using Ensembl or GENCODE.
  # Input: ORF table and optional GENCODE UniProt mapping file.
  # Output: ORF table with UniProt protein IDs linked to gene ID.
  
  annotation_provider <- tolower(annotation_provider)  # Standardize input case
  
  if (annotation_provider == "ensembl") {
    # Connect to Ensembl biomart database
    mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Retrieve Ensembl Gene ID to UniProt ID mapping
    annotLookup <- biomaRt::getBM(
      attributes = c('ensembl_gene_id', 'uniprot_gn_id'),
      mart = mart
    ) %>%
      dplyr::group_by(ensembl_gene_id) %>%
      dplyr::mutate(uniprot_gn_ids = paste0(uniprot_gn_id, collapse = ";")) %>%  # Combine multiple UniProt IDs per gene
      dplyr::ungroup() %>%
      dplyr::select(ensembl_gene_id, uniprot_gn_ids) %>%
      dplyr::distinct()
    
    # Merge with ORF table
    orf_table <- orf_table %>%
      dplyr::left_join(annotLookup, by = c("gene_id" = "ensembl_gene_id"))
    
  } else if (annotation_provider == "gencode") {
    # Use user-provided GENCODE mapping file
    gencode_conversion <- read.delim(gencode_uniprot_file, header = FALSE) %>%
      dplyr::select(-V3) %>%  # Remove unused column
      dplyr::rename(transcript_id = V1, uniprot_gn_id = V2) %>%
      
      # Map transcript IDs to gene IDs in ORF table
      dplyr::left_join(orf_table %>%
                         dplyr::select(gene_id, transcript_id) %>%
                         dplyr::distinct(), by = "transcript_id") %>%
      
      # Ensure each gene has at least one UniProt entry
      dplyr::filter(complete.cases(.)) %>%
      dplyr::select(-transcript_id) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(uniprot_gn_ids = paste0(uniprot_gn_id, collapse = ";")) %>%  # Combine UniProt IDs per gene
      dplyr::ungroup() %>%
      dplyr::select(gene_id, uniprot_gn_ids) %>%
      dplyr::distinct()
  }
}

calculate_protein_similarity <- function(i) {
  # Calculate protein similarity between uniprot sequence and protein sequence of
  # annotated ORF
  # Input: single ORF protein sequence
  #        uniprot protein sequence
  # Output: similarity score between ORF and uniprot protein
  
  # Get UniProt IDs associated with the ORF
  prot_ids <- orf_table$uniprot_gn_ids[i]
  if (is.na(prot_ids)) return(NA)
  
  orf_sequence <- Biostrings::AAString(orf_table$Protein[i])
  prot_ids_split <- strsplit(prot_ids, ";")[[1]]
  prot_ids_split <- prot_ids_split[prot_ids_split %in% names(uniprot_fasta)]
  fasta_entries <- uniprot_fasta[prot_ids_split]
  fasta_entries <- fasta_entries[!sapply(fasta_entries, is.null)]

  if (length(fasta_entries) == 0) return(NA)

  alignment <- pwalign::pairwiseAlignment(pattern = fasta_entries, subject = orf_sequence, type = "local")
  similarity_score <- max(nmatch(alignment) / width(fasta_entries)) * 100
  return(similarity_score)
}


# Load correct ORF information depending on used ORFcaller
if ( tolower(orfcaller) == "orfquant") {
  
  prep_orfs <- prepare_orfquant(orfquant_orfs_loc = orfs_loc)
  
} else if ( tolower(orfcaller) == "price") {
  
  prep_orfs <- prepare_price(price_orfs_loc = orfs_loc,
                             tx2gene = tx2gene)
}

# Check annotation style
restyled_orfs <- check_annot_style(prep_orfs[[1]],
                                   annotated_gen = cds_gene)

# Find new ORF and annotated ORF location similarities
orf_cds_sim <- check_orf_cds_similarity(orf_ranges = restyled_orfs,
                                        orf_table = prep_orfs[[2]],
                                        annotated_gen = cds_gene,
                                        annotated_gen_unlist = cds_gene_unlist)

# Annotate ORF table with new found annotations
orf_table <- annotate_new_orfs(orf_ranges = restyled_orfs,
                               orf_table = prep_orfs[[2]],
                               cds_matches_grl = orf_cds_sim,
                               orf_caller = orfcaller)

# Annotate ORFs with uniprot protein IDs using either gencode
# conversion table or Ensembl biomart
orf_table <- annotate_uniprot_id(annotation_provider,
                                 gencode_uniprot_file,
                                 orf_table)

# Load UniProt Reference Proteome fasta files
uniprot_fasta <- rtracklayer::import(uniprot_protein_fasta_loc,
                                     type = "AA")

# Extract uniprot names from fasta headers
# Might have to be changed when the input file is different
names(uniprot_fasta) <- sapply(names(uniprot_fasta), function(x) {
  strsplit(x, "\\|")[[1]][2]
})


# Calculate protein similarity in parallel
orf_table$similarity_score <- parallel::mclapply(1:nrow(orf_table), calculate_protein_similarity, mc.cores = cpus)
orf_table$similarity_score <- as.numeric(orf_table$similarity_score)

# Write resulting table to disk
write.table(orf_table, file = paste("/hpc/pmc_vanheesch/projects/evanderwerf/nextflow_pipelines/riboseq_pipeline/annotation_of_expression_data",paste(orfcaller,"orfs.csv", sep = "_"),sep = "/"),
            quote = F,
            row.names = F,
            sep = ",")
