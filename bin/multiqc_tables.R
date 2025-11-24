#!/usr/bin/env Rscript

### Load libraries
suppressPackageStartupMessages({
    library(stringr)
    library(ggplot2)
    library(ggsci)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(RColorBrewer)
    library(RiboseQC)
})

### Catch arguments
arguments <- commandArgs(trailingOnly = TRUE)

# Get first argument to decide correct function to run
function_to_run <- arguments[1] 

### Load functions
multiqc_canoncial_count_table <- function(arguments) {

    orf_table_csv = arguments[2]
    orf_expression_csv = arguments[3]
    
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
        dplyr::count(variable, orf_category_group, name = "count") %>%
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

multiqc_riboseqc_tables <- function(arguments){
    # Load RiboseQC output files and initialize color palette
    input_files <- arguments[-1]

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
        dplyr::filter(read_length %in% 28:31, comp == "nucl") %>%
        dplyr::select(sample_id, read_length, frame_preference) %>%
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
        dplyr::select(sample_id, all_of(cols_sorted))
    
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



# Periodicity plots load input data
load_psite_data <- function(riboseqc_files, read_length = "29") {
    # Define CDS region (based on your metagene plot positions)
    cds_start <- 51  # position of start codon
    cds_stop <- 150  # position of stop codon
    
    combined_signal <- numeric(199)
    sample_frames <- list()
    
    for(filepath in riboseqc_files) {
        temp_env <- new.env()
        load(filepath, envir = temp_env)
        
        psite_data <- temp_env$res_all[["profiles_P_sites"]][["P_sites_subcodon"]][["nucl"]][[read_length]]
        signal <- colSums(as.matrix(psite_data))
        
        # Add to combined signal for metagene
        combined_signal <- combined_signal + signal
        
        # Calculate frame percentages only for CDS region
        cds_signal <- signal[cds_start:cds_stop]
        frame_counts <- sapply(1:3, function(i) {
        frame_indices <- which(rep(c(2,3,1), length.out = length(cds_signal)) == i)
        sum(cds_signal[frame_indices])
        })
        sample_frames[[basename(filepath)]] <- frame_counts / sum(frame_counts) * 100
        
        rm(temp_env)
        gc()
    }
    
    return(list(
        metagene = combined_signal,
        periodicity = do.call(rbind, sample_frames)
    ))
}

# Plotting function for metagene
plot_psite_metagene <- function(signal_data, title = "Combined P-site metagene profile (29nt reads)") {
    plot_df <- data.frame(
        pos = seq_along(signal_data),
        value = signal_data,
        frame = factor(paste("Frame", rep(c(2,3,1), length.out = length(signal_data))),
                levels = paste("Frame", 1:3))
    )
    metagene_plot <- ggplot(plot_df, aes(x = pos, y = value, fill = frame, color = frame)) +
        geom_bar(stat = "identity", width = 0.5, position = "identity") +  # width=1 makes bars touch
        scale_fill_manual(values = c("#985143", "#32936F", "#1B365D")) +
        scale_color_manual(values = c("#985143", "#32936F", "#1B365D")) +
        scale_x_continuous(
            breaks = c(1, 26, 51, 84, 117, 150, 176, 200),
            labels = c("TSS", "", "start\ncodon", "", "", "stop\ncodon", "", "TES"),
            expand = c(0, 0)  # removes spacing at edges
        ) +
        scale_y_continuous(expand = c(0, 0)) +  # removes spacing at bottom
        labs(
            x = "Position (nucleotide resolution)",
            y = "P-site count",
            title = title
        ) +
        theme_minimal() +
        theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 16, face = "bold"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()
        )
    ggsave(filename = "Metagene_profile_combined_mqc.png",
        plot = metagene_plot,
        width = 7,
        height = 5,
        dpi = 150)
    metagene_plot

}

# Plotting function for periodicity
plot_periodicity <- function(periodicity_data) {
    # Convert to long format for plotting
    colnames(periodicity_data) <- c("Frame 3", "Frame 1", "Frame 2")
    plot_df <- as.data.frame(periodicity_data) %>%
        tibble::rownames_to_column("sample") %>%
        tidyr::pivot_longer(-sample, 
                            names_to = "frame", 
                            values_to = "percentage") %>%
        mutate(frame = factor(frame))
    
    # Calculate summary statistics
    summary_stats <- plot_df %>%
        group_by(frame) %>%
        summarise(
        mean_pct = mean(percentage),
        sd_pct = sd(percentage)
        )
    
    periodicity_plot <- ggplot() +
        geom_jitter(data = plot_df, 
                    aes(x = frame, y = percentage, color = frame),
                    width = 0.2, 
                    alpha = 1,
                    size = 2) +
        geom_bar(data = summary_stats,
                aes(x = frame, y = mean_pct, fill = frame),
                stat = "identity",
                alpha = 0.5,
                width = 0.9) +
        geom_errorbar(data = summary_stats,
                    aes(x = frame, 
                        ymin = mean_pct - sd_pct,
                        ymax = mean_pct + sd_pct),
                    width = 0.2) +
        scale_y_continuous(
        limits = c(0, max(summary_stats$mean_pct + summary_stats$sd_pct) * 1.1),
        expand = expansion(mult = c(0, 0)),
        breaks = seq(0, 100, by = 10)  # for every 10%
        ) +
        scale_fill_manual(values = c("#985143", "#32936F", "#1B365D")) +
        scale_color_manual(values = c("#985143", "#32936F", "#1B365D")) +
        labs(
        x = "Reading frame",
        y = "Ribosome P-sites (%)"
        ) +
        theme_classic() +
        theme(
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()
        )
    
    ggsave(filename = "Periodicity_bar_combined_mqc.png",
            plot = periodicity_plot,
            width = 7,
            height = 5,
            dpi = 110)
    periodicity_plot
}


# Run specific function by looking at first args value
if (function_to_run == "riboseqc_tables"){
    multiqc_riboseqc_tables(arguments)

} else if(function_to_run == "canonical_count"){
    multiqc_canoncial_count_table(arguments)

} else if(function_to_run == "periodicity_plot"){

    riboseqc_files <- as.character(arguments[-1])
    # Load riboseqc data  
    psite_data <- load_psite_data(riboseqc_files)
    # Create plots as needed
    metagene_plot <- plot_psite_metagene(psite_data$metagene)
    periodicity_plot <- plot_periodicity(psite_data$periodicity)

} else {
    print(paste0("No function found that matches: ", function_to_run))
}



