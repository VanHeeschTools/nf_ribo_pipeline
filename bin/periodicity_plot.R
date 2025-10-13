#!/usr/bin/env Rscript
suppressPackageStartupMessages({
    library(RiboseQC)
    library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
riboseqc_files <- files <- as.character(args)

# Main function to load and process data
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
    ggplot(plot_df, aes(x = pos, y = value, fill = frame, color = frame)) +
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

    ggplot() +
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
}

# Run function

# Load data once
psite_data <- load_psite_data(riboseqc_files)

# Create plots as needed
metagene_plot <- plot_psite_metagene(psite_data$metagene)
ggsave(filename = "Metagene_profile_combined_mqc.png",
    width = 7,
    height = 5,
    dpi = 150)
metagene_plot

periodicity_plot <- plot_periodicity(psite_data$periodicity)
ggsave(filename = "Periodicity_bar_combined_mqc.png",
    width = 7,
    height = 5,
    dpi = 110)
periodicity_plot
