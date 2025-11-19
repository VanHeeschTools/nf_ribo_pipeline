multiqc_canoncial_count_table <- function(orf_table_csv, orf_expression) {
    # Load libraries
    suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(stringr)
        library(RColorBrewer)
    })

    # Load intput csv
    orf_table <- read.csv(orf_table_csv)
    orf_expression <- read.csv(orf_expression_csv)
    
    # Join and clean expression column names
    expression_cols <- setdiff(names(orf_expression), "orf_id")
    cleaned_cols <- str_remove(expression_cols, "^X")
    
    orf_table <- orf_table %>%
        left_join(orf_expression, by = "orf_id") %>%
        rename_with(~ cleaned_cols, all_of(expression_cols))
    
    # Convert expression columns to logical TRUE/FALSE (PPM >= 1)
    orf_table <- orf_table %>%
        mutate(across(all_of(cleaned_cols),
                    ~ as.numeric(.) >= 1),
            orf_category_group = if_else(
                orf_biotype_single %in% c("ORF-annotated", "NC-variant"),
                "Canonical", "Non-canonical"
            ))
    
    # Count TRUE per sample and category
    df_counts_wide <- orf_table %>%
        pivot_longer(
        cols = all_of(cleaned_cols),
        names_to = "variable",
        values_to = "value"
        ) %>%
        filter(value) %>%
        count(variable, orf_category_group, name = "count") %>%
        pivot_wider(
        names_from = orf_category_group,
        values_from = count,
        values_fill = 0
        )
    
    # Write output to txt file
    write.table(
        df_counts_wide,
        "canonical_orf_counts_mqc.txt",
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


multiqc_riboseqc_tables <- function(input_files){
    # Load libraries
    suppressPackageStartupMessages({
        library(stringr)
        library(ggplot2)
        library(ggsci)
        library(dplyr)
        library(tidyr)
        library(scales)
    })


    # Load RiboseQC output files and initialize color palette
    input_files <- commandArgs(trailingOnly = TRUE)
    
    if (length(input_files) == 0) stop("No input files provided.")
    
    # Initialize empty data frames for collecting statistics from each sample
    summary_P_sites_df <- data.frame()
    summary_reads_df <- data.frame()
    inframe_df <- data.frame()
    read_cats_df <- data.frame()
    cds_reads_df <- data.frame()
    
    # Process each RiboseQC results file
    for (fname in input_files) {
        parts <- str_split(basename(fname), "_")[[1]]
        sample_id <- paste(parts[1:(length(parts) - 3)], collapse = "_")
        
        message("Loading ", sample_id)
        load(fname)
        
        # Extract data from loaded RiboseQC result
        summary_P_sites_sample <- as.data.frame(res_all$summary_P_sites)
        summary_P_sites_sample$sample_id <- sample_id
        
        summary_reads_sample <- as.data.frame(t(colSums(as.data.frame(res_all$read_stats$reads_summary_unq$nucl))))
        summary_reads_sample$sample_id <- sample_id
        
        inframe_sample <- as.data.frame(t(res_all$selection_cutoffs$analysis_frame_cutoff$nucl$all$frames_res))
        rownames(inframe_sample) <- sample_id
        
        read_cats_sample <- as.data.frame(t(rowSums(as.data.frame(res_all$read_stats$reads_summary$nucl))))
        rownames(read_cats_sample) <- sample_id
        
        cds_reads_sample <- data.frame(reads = sum(res_all$read_stats$counts_cds_genes_unq$reads))
        rownames(cds_reads_sample) <- sample_id
        
        # Merge into master data frames
        summary_P_sites_df <- rbind(summary_P_sites_df, summary_P_sites_sample)
        summary_reads_df <- bind_rows(summary_reads_df, summary_reads_sample)
        inframe_df <- rbind(inframe_df, inframe_sample)
        read_cats_df <- rbind(read_cats_df, read_cats_sample)
        cds_reads_df <- rbind(cds_reads_df, cds_reads_sample)
    }
    
    # Frame Preference Table
    summary_P_sites_df_s <- summary_P_sites_df  
    
    # Frame preference
    df_wide <- summary_P_sites_df_s %>%
        filter(read_length %in% 28:31, comp == "nucl") %>%
        select(sample_id, read_length, frame_preference) %>%
        pivot_wider(
        id_cols = sample_id,
        names_from = read_length,
        values_from = frame_preference,
        names_prefix = "percentage_",
        values_fill = NA
        )
    # Sort the "percentage_" columns numerically
    cols <- setdiff(names(df_wide), "sample_id")
    cols_sorted <- cols[order(as.numeric(gsub("\\D", "", cols)))]
    
    # Reorder and show result
    df_wide <- df_wide %>%
        select(sample_id, all_of(cols_sorted))
    
    write.table(
        df_wide,
        file = "inframe_percentages_mqc.txt",  
        sep = "\t",                        
        quote = FALSE,                    
        row.names = FALSE           
    )
    
    # Obtain counts of every read category for every input sample
    read_cats_raw <- read_cats_df
    read_cats_raw$Sample <- rownames(read_cats_raw)
    read_cats_raw <- read_cats_raw[, c(ncol(read_cats_raw), 1:(ncol(read_cats_raw) - 1))]
    
    header <- paste(colnames(read_cats_raw), collapse = "\t")
    lines <- apply(read_cats_raw, 1, function(row) paste(row, collapse = "\t"))
    
    writeLines(c(
        header,
        lines
    ), con = "riboseqc_read_categories_counts_mqc.txt")
}


