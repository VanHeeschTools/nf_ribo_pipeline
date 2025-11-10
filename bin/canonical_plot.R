#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(tidyr)
    library(stringr)
    library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
orf_table_csv <- args[1]
orf_expression_csv <- args[2]

orf_table <- read.csv(orf_table_csv)
orf_expression <- read.csv(orf_expression_csv)

# Join expression data to ORF table
orf_table <- orf_table %>% 
    left_join(orf_expression, by = "orf_id")

# -------------------------------
# Define canonical ORF categories
# -------------------------------
annotated_vector <- c("ORF-annotated", "NC-variant")

# -------------------------------
# Select all expression columns added by orf_expression (except orf_id)
# -------------------------------
expression_cols <- setdiff(names(orf_expression), "orf_id")

# -------------------------------
# Clean column names (remove X and keep before first "_")
# -------------------------------
clean_names <- function(x) {
    x %>%
        str_remove("^X") # remove leading X
}

names(orf_table)[match(expression_cols, names(orf_table))] <- clean_names(expression_cols)
expression_cols <- clean_names(expression_cols)  # update vector

# -------------------------------
# Convert expression columns to logical TRUE/FALSE (PPM >= 1)
# -------------------------------
orf_table <- orf_table %>%
    mutate(across(all_of(expression_cols),
                ~ if_else(!is.na(suppressWarnings(as.numeric(.))),
                            as.numeric(.) >= 1,   # TRUE if PPM >= 1
                            FALSE))) %>%
    # Define ORF category group: Canonical vs Non-Canonical
    mutate(orf_category_group = if_else(orf_biotype_single %in% annotated_vector,
                                    "Canonical", "Non-canonical"))

# -------------------------------
# Reshape expression columns to long format
# -------------------------------
orf_table_long <- orf_table %>%
    pivot_longer(
        cols = all_of(expression_cols),
        names_to = "variable",
        values_to = "value"
    )

# -------------------------------
# Count TRUE values per variable grouped by ORF category
# -------------------------------
orf_table_counts <- orf_table_long %>%
    filter(value) %>%   # keep only TRUEs
    group_by(variable, orf_category_group) %>%
    summarise(count = n(), .groups = "drop")

df_counts_wide <- orf_table_counts %>%
    tidyr::pivot_wider(
        names_from = orf_category_group,
        values_from = count,
        values_fill = 0
)

write.table(
    df_counts_wide,
    file = "canonical_orf_counts_mqc.txt",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
