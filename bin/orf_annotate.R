suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(GenomicRanges)
  library(rtracklayer)
  library(stringr)
  library(AnnotationDbi)
  library(Biostrings)
  library(parallel)
})

# This script creates new annotations for PRICE and ORFquant ORFs based on the annotated
# reference ORFs

args <- commandArgs(trailingOnly = TRUE)
orfs_loc <- args[1]
txdb_loc <- args[2]
orfcaller <- args[3]
annotation_provider <- args[4]
gencode_uniprot_file <- args[5]
uniprot_protein_fasta_loc <- args[6]
cpus <- args[7]

# Extract CDS regions of annotated genes from txdb
txdb <- AnnotationDbi::loadDb(txdb_loc)
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
cds_gene <- GenomicFeatures::cdsBy(txdb, "gene")
cds_gene_unlist <- unlist(cds_gene)

# FUNCTIONS
prepare_price <- function(price_orfs_loc,
                          tx2gene) {

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
    dplyr::select(!(gene_id))

    return(list(orf_ranges,
                price_orf_df))
}

prepare_orfquant <- function(orfquant_orfs_loc) {

  # This functions imports the ORFquant object and outputs
  # a list of ORF ranges (orf_ranges), 
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

check_annot_style <- function(orf_ranges,
                              annotated_gen) {
  
  # Convert the seqlevels style of a GRange so that they can be compared with one
  # another. So far only checks for NCBI and Ensembl.
  # Input: orf_ranges and CDS definitions of the ORF, only checks the seqlevels of
  #        these GRanges
  # Output: orf_ranges with updated seqlevels
  
  # This is quite hard-coded, I would recommend checking the seqlevels of the
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

check_orf_cds_similarity <- function(orf_ranges,
                                     orf_table,
                                     annotated_gen,
                                     annotated_gen_unlist) {

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

annotate_new_orfs <- function(orf_ranges,
                              orf_table,
                              cds_matches_grl,
                              orf_caller) {

  # Compare CDS vs new ORF, and make a new annotation based on a virtual longest
  # possible ORF of the annotated set.
  # Input: orf_ranges (Granges of orf_caller ORFs)
  #        orf_table (metadata of orf_caller ORFs)
  #        cds_matches_grl ()
  #        orf_caller (This is important as ORFquant does not include the STOP codon in the sequence)
  # Output: orf_table with new columns showing start codon, orf coord similarity and the new
  #         orf annotation type.

  # ORFquant does not include stops in the ORF, so we trim the CDS matches

  if (tolower(orf_caller) == "orfquant") {
    cds_range_similarity <- width(range(orf_ranges)) / (width(range(cds_matches_grl))-3)
    cds_strand <- ifelse(elementNROWS(cds_matches_grl) > 0,
                         as.character(unique(strand(cds_matches_grl))),
                         NA)

    # Extract strand values from S4 object
    orf_strand <- as.character(unlist(runValue(strand(orf_ranges))))
    orf_start <- ifelse(orf_strand == "+", 
                        min(start(orf_ranges)),
                        max(end(orf_ranges)))

    # pay attention to see whether your ORF caller includes the stop or not
    orf_stop <- ifelse(orf_strand == "+",
                        max(end(orf_ranges)),
                        min(start(orf_ranges)))

    ann_start <- ifelse(cds_strand == "+",
                        min(start(cds_matches_grl)),
                        max(end(cds_matches_grl)))
    # pay attention to see whether your ORF caller includes the stop or not
    ann_stop <- ifelse(cds_strand == "+",
                        max(end(cds_matches_grl)) - 3,
                        min(start(cds_matches_grl)) + 3)

    orf_category <- rep("Unknown", length(orf_ranges))

  } else if (tolower(orf_caller) == "price") {

    cds_range_similarity <- width(range(orf_ranges)) / (width(range(cds_matches_grl)))
    cds_strand <- ifelse(elementNROWS(cds_matches_grl) > 0,
                         as.character(unique(strand(cds_matches_grl))),
                         NA)

    # Extract strand values from S4 object
    orf_strand <- as.character(unlist(runValue(strand(orf_ranges))))
    orf_start <- ifelse(orf_strand == "+", 
                        min(start(orf_ranges)),
                        max(end(orf_ranges)))

    # pay attention to see whether your ORF caller includes the stop or not
    orf_stop <- ifelse(orf_strand == "+",
                        max(end(orf_ranges)),
                        min(start(orf_ranges)))

    ann_start <- ifelse(cds_strand == "+",
                        min(start(cds_matches_grl)),
                        max(end(cds_matches_grl)))
    # pay attention to see whether your ORF caller includes the stop or not
    ann_stop <- ifelse(cds_strand == "+",
                        max(end(cds_matches_grl)),
                        min(start(cds_matches_grl)))

    orf_category <- rep("Unknown", length(orf_ranges))
  }

  # Positive strand
  pos_strand_idx <- orf_strand == "+"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start == ann_start] <- "ORF_annotated"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start < ann_start] <- "N_extension"
  orf_category[pos_strand_idx & orf_stop == ann_stop & orf_start > ann_start] <- "N_truncation"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop < ann_stop] <- "overl_uORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop < ann_start] <- "uORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop > ann_stop] <- "NC_extension"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop > ann_stop] <- "overl_dORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_stop & orf_stop > ann_stop] <- "dORF"
  orf_category[pos_strand_idx & orf_stop != ann_stop & orf_start > ann_start & orf_stop < ann_stop] <- "nested_ORF"
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
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start < ann_start & orf_stop > ann_stop] <- "nested_ORF"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop > ann_stop] <- "C_truncation"
  orf_category[neg_strand_idx & orf_stop != ann_stop & orf_start == ann_start & orf_stop < ann_stop] <- "C_extension"

  # The previous code is for CDS matches, we do not have the info for novel ORFs
  orf_category[lengths(cds_matches_grl) == 0] <- "novel"
  cds_range_similarity[lengths(cds_matches_grl) == 0] <- NA

  if (tolower(orf_caller == "orfquant")) {
    new_category_df <- data.frame(orf_id = names(orf_ranges),
                                    start_dif = abs(orf_start - ann_start),
                                    cds_range_similarity = as.numeric(cds_range_similarity),
                                    stop_same = ifelse(orf_stop == ann_stop,
                                    TRUE, FALSE),
                                    orf_category_new = orf_category,
                                    start_codon = "ATG",
    # range similarity can be bigger due to new ORF being bigger than the annotated ORF
                                    orf_cds = ifelse(cds_range_similarity >= 0.9 &
                                                    cds_range_similarity <= 1.1,
                                                    TRUE, FALSE)) %>%
        dplyr::mutate(start_check = ifelse(start_dif > 0,
                                        TRUE, FALSE),
    # Simple checks for similarity
                    same_as_cds = ifelse(start_dif < 99 &
                                        start_check == TRUE &
                                        stop_same == TRUE &
                                        orf_cds == TRUE,
                                        TRUE, FALSE))

    # Add everything to a complete table
    orf_table <- orf_table %>%
    dplyr::left_join(new_category_df, by = c("ORF_id_tr" = "orf_id"))

    return(orf_table)

  } else if (tolower(orf_caller) == "price") {
    new_category_df <- data.frame(orf_id = names(orf_ranges),
                                  start_dif = abs(orf_start - ann_start),
                                  cds_range_similarity = as.numeric(cds_range_similarity),
                                  stop_same = ifelse(orf_stop == ann_stop,
                                  TRUE, FALSE),
                                  orf_category_new = orf_category,
                                  start_codon = stringr::str_split_i(names(orf_ranges),"__", i = 5),
    # range similarity can be bigger due to new ORF being bigger than the annotated ORF
                                  orf_cds = ifelse(cds_range_similarity >= 0.9 &
                                                   cds_range_similarity <= 1.1,
                                                   TRUE, FALSE)) %>%
    dplyr::mutate(start_check = ifelse(start_dif > 0,
                                       TRUE, FALSE),
    # Simple checks for similarity
                  same_as_cds = ifelse(start_dif < 99 &
                                       start_check == TRUE &
                                       stop_same == TRUE &
                                       orf_cds == TRUE &
                                       !(start_codon == "ATG"),
                                       TRUE, FALSE))

    # Add everything to a complete table
    orf_table <- orf_table %>%
    dplyr::select(!start_codon) %>%
    dplyr::left_join(new_category_df, by = c("name" = "orf_id"))

    return(orf_table)
  }
}

annotate_uniprot_id <- function(orf_table) {

  # Annotate ORFs with uniprot ID where possible using biomart
  # Input: ORF table
  # Output: ORF table with uniprot protein IDs linked to gene ID

  if ( tolower(annotation_provider) == "ensembl") {

  mart <- biomaRt::useEnsembl("ENSEMBL_MART_ENSEMBL")
  mart <- biomaRt::useDataset('hsapiens_gene_ensembl', mart)
  annotLookup <- biomaRt::getBM(
  mart <- mart,
  attributes = c(
    'ensembl_gene_id',
    'uniprot_gn_id'),
    uniqueRows = FALSE) %>%
  dplyr::group_by(ensembl_gene_id) %>%
  # Add all possible uniprot IDs per gene together in single value
  dplyr::mutate(uniprot_gn_ids = paste0(uniprot_gn_id, collapse = ";")) %>%
  ungroup() %>%
  dplyr::select(c(ensembl_gene_id, uniprot_gn_ids)) %>%
  dplyr::distinct()

  orf_table <- orf_table %>%
  dplyr::left_join(annotLookup, by = c("gene_id" = "ensembl_gene_id"))

  } else if (tolower(annotation_provider == "gencode")) {
    
    # Use file from gencode provided by the user
    # Have to link transcript IDs to gene IDs in the ORF table
    gencode_conversion <- read.delim(gencode_uniprot_file, header = F) %>%
      dplyr::select(-V3) %>%
      dplyr::rename(transcript_id = V1,
                    uniprot_gn_id = V2) %>%
      dplyr::left_join(orf_table %>% 
                         dplyr::select(gene_id,transcript_id) %>% 
                         dplyr::distinct(), by = "transcript_id") %>%
      # Not all genes are present in the called ORFs. This will make sure that
      # we have at least 1 uniprot per gene of each called ORF from the software
      dplyr::filter(complete.cases(.)) %>%
      dplyr::select(-transcript_id) %>%
      dplyr::distinct() %>%
      dplyr::group_by(gene_id) %>%
      # Add all possible uniprot IDs per gene together in single value
      dplyr::mutate(uniprot_gn_ids = paste0(uniprot_gn_id, collapse = ";")) %>%
      dplyr::ungroup() %>%
      dplyr::select(c(gene_id, uniprot_gn_ids)) %>%
      dplyr::distinct()
  }
}

# Calculate protein similarity between uniprot sequence and protein sequence of
# annotated ORF
# Input: single ORF protein sequence
#        uniprot protein sequence
# Output: similarity score between ORF and uniprot protein
process_sequence <- function(i) {
  prot_ids <- orf_table$uniprot_gn_ids[i]
  if (is.na(prot_ids)) return(NA)
  
  orf_sequence <- Biostrings::AAString(orf_table$Protein[i])
  prot_ids_split <- strsplit(prot_ids, ";")[[1]]
  
  prot_ids_split <- prot_ids_split[prot_ids_split %in% names(uniprot_fasta_names)]
  fasta_entries <- uniprot_fasta_names[prot_ids_split]
  fasta_entries <- fasta_entries[!sapply(fasta_entries, is.null)]
  
  if (length(fasta_entries) == 0) return(NA)
  
  # prot_seqs <- AAStringSet(unlist(fasta_entries))
  alignment <- Biostrings::pairwiseAlignment(pattern = fasta_entries, subject = orf_sequence, type = "local")
  similarity_score <- max(nmatch(alignment) / width(fasta_entries)) * 100
  
  return(as.numeric(similarity_score))
}

# Load ORF information
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

# Annotate ORF table with newfound annotations
orf_table <- annotate_new_orfs(orf_ranges = restyled_orfs,
                               orf_table = prep_orfs[[2]],
                  cds_matches_grl = orf_cds_sim,
                  orf_caller = orfcaller)

# Annotate ORFs with uniprot protein IDs using either gencode
# conversion table or Ensembl biomart
orf_table <- annotate_uniprot_id("gencode",
                                 gencode_uniprot_file,
                                 orf_table)

# Load UniProt Reference Proteome fasta files
uniprot_fasta <- rtracklayer::import(uniprot_fasta_loc,
                                     type = "AA")

# Change the names to more easily digestible by R
# Might have to be changed when the input file is different
names(uniprot_fasta) <- sapply(names(uniprot_fasta), function(x) {
  strsplit(x, "\\|")[[1]][2]
})

# Do the calculations in parallel using MCLAPPLY
orf_table <- orf_table %>%
  dplyr::mutate(similarity_score = parallel::mclapply(1:nrow(orf_table),
                                                      process_sequence,
                                                      mc.cores = cpus))

# Write resulting table to disk
write.table(orf_table, file = paste(orf_caller,"orfs.csv", sep = "_"),
            quote = F,
            row.names = F,
            sep = ",")
