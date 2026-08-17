#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(rtracklayer)
})

# Load input data
args <- commandArgs(trailingOnly = TRUE)
orfcaller_tables <- args[1]

#' Turn the harmonised ORF table into a gtf-like file
#' 
#' @param sorted_df data.frame produced by orf_filter()
convert_to_gtf <- function(sorted_df, gtf_output_file) {

    orf_table <- read.csv(sorted_df)

    # Find all columns that show if ORF is found in caller
    orf_unique <- sorted_df %>%
        dplyr::select(orf_id, dplyr::all_of(found_cols)) %>%
        dplyr::distinct(orf_id, .keep_all = TRUE)

    orf_found_attrs <- data.frame(
        orf_id = orf_unique$orf_id,
        found_attrs = do.call(
            paste0,
            base::lapply(found_cols, function(col) paste0(col, ' "', orf_unique[[col]], '"; '))
        ),
        stringsAsFactors = FALSE
    )

    # Handle transcript rows
    message("Building transcript-level rows")
    transcripts <- orf_table %>%
    left_join(orf_found_attrs, by = "orf_id") %>%
    mutate(
        source = "Harmonised",
        start = as.integer(orf_start), 
        end   = as.integer(orf_end),
        feature    = "transcript",
        score      = ".",
        frame      = ".",
        attributes = paste0('transcript_id "', orf_id, '"; gene_id "', 
                            gene_id, '"; gene_name "', gene_name, 
                            '"; gene_biotype "', gene_biotype,
                            '"; ORF_id "', orf_id,
                            '"; ORF_biotype_single "', orf_biotype_single, '"; ',
                            found_attrs)
    ) %>%
    # Orf_id will be used to join with cds rows, and is removed afterwards
    dplyr::select(chr, source, feature, start, end, 
                    score, strand, frame, attributes, orf_id) %>%
    arrange(chr, start) # Sort based on genomic location

    message("  -> Built ", nrow(transcripts), " transcript rows")

    message("Building CDS rows (splitting multi-exon start/end fields)")
    cds <- orf_table %>%
    mutate(
        source = "Harmonised",
        starts = strsplit(as.character(starts), "_"),
        ends   = strsplit(as.character(ends), "_")
    ) %>%
    unnest(c(starts, ends)) %>%
    left_join(orf_found_attrs, by = "orf_id") %>% 
    mutate(
        start = as.integer(starts),
        end   = as.integer(ends)
    ) %>%
    mutate(
        feature    = "CDS",
        score      = ".",
        frame      = ".",
        attributes = paste0('transcript_id "', orf_id, '"; gene_id "', 
                            gene_id, '"; gene_name "', gene_name, 
                            '"; gene_biotype "', gene_biotype,
                            '"; ORF_id "', orf_id,
                            '"; ORF_biotype_single "', orf_biotype_single, '"; ',
                            found_attrs)
    ) %>%
    dplyr::select(chr, source, feature, start, end, 
                    score, strand, frame, attributes, orf_id)

    message("  -> Built ", nrow(cds), " CDS rows")
    message("Merging transcript and CDS rows per ORF (", nrow(transcripts), " transcripts), this might take a while")

    # Combine transcript rows with their corresponding CDS rows
    # The CDS rows are sorted by strand
    n_tx <- nrow(transcripts) 
    cds_by_orf <- split(cds, cds$orf_id)

    gtf_list <- lapply(1:n_tx, function(i) {
        tx <- transcripts[i, ]
        tx_cds <- cds_by_orf[[tx$orf_id]]
        if (is.null(tx_cds)) tx else bind_rows(tx, tx_cds)
    })

    message("  -> Finished merging all transcripts")
    message("Combining all rows into final GTF table...")

    # Combine all transcript+CDS groups into a single tibble
    gtf_out <- do.call(rbind, gtf_list) %>%
    dplyr::select(-orf_id) %>% # remove helper column used for grouping
    mutate(chr = as.character(chr)) # convert factor back to character

    message("  -> Final GTF has ", nrow(gtf_out), " total rows")
    message("Writing output GTF to: ", gtf_output_file)

    # Write output gtf file
    write.table(gtf_out, file = gtf_output_file, sep = "\t", quote = FALSE, 
                col.names = FALSE, row.names = FALSE)
    message("Finished conversion script")
}


# Run conversion script
convert_to_gtf(orfcaller_tables, "harmonised_orf_table.gtf")