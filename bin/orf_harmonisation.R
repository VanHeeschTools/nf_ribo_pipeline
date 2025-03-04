#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(rtracklayer)
library(Biostrings)
library(ggplot2)
paths <- c("bsgenome_install/", .libPaths())
.libPaths(paths)

library(BSgenome.Homo.sapiens)
genome <- BSgenome.Homo.sapiens


basedir <- ""
txdb_loc <- ""
orfquant_table <- ""
price_table <- ""
orfquant_ppm <- ""
price_ppm <- ""


orfquant_ppm_table <- read.csv(orfquant_ppm, header = TRUE,  row.names = NULL)
orfquant_ppm_table$number_patient_samples <- apply(orfquant_ppm_table[, -1], 1, function(x) sum(x >= 1))

price_ppm_table <- read.csv(price_ppm, header = TRUE,  row.names = NULL)
price_ppm_table$number_patient_samples <- apply(price_ppm_table[, -1], 1, function(x) sum(x >= 1))

# Also calculate ppm mean and ppm median?

min_number_of_samples <- 2
min_similarity_score <- 80



# TRUE DATA LOAD ---------

# Load ORF annotation table from PRICE and ORFquant
# TODO: Need to double check for input consistency between old and new code
price_orfs <- read.delim(price_table, sep = ",") %>%
  # Join the number of patient samples from the price_ppm_table
  dplyr::left_join(
    price_ppm_table %>%
      dplyr::select(orf_id, number_patient_samples),
    by = c("name" = "orf_id")
  ) %>%
  
  # Create a new column 'ORF_ranges' that combines genomic coordinates
  # and mark ORFs as translated if detected in more than the minimum number of samples
  dplyr::mutate(
    ORF_ranges = paste0(seqnames, ":", start, "-", end),
    translated = ifelse(number_patient_samples > min_number_of_samples, TRUE, FALSE)
  ) %>%
  
  # Remove unnecessary columns
  dplyr::select(-c("start_check","orf_cds","start_dif","stop_same","seqnames", "start", "end", "width", "strand")) %>%
  
  # Rename column 'name' to 'orf_id' for consistency
  dplyr::rename(orf_id = name) %>%
  
  # Add a source column for better tracking
  dplyr::mutate(source = "PRICE") %>%
  mutate(P_sites_raw = NA, 
         P_sites_raw_uniq = NA, 
         ORF_category_Tx = NA,
         ORF_category_Gen = NA
  )
  
  
orfquant_orfs <- read.delim(orfquant_table, sep = ",") %>%
  
  # Add the number of samples with a PPM higher than 1
  dplyr::left_join(
    orfquant_ppm_table %>%
      dplyr::select(orf_id, number_patient_samples),
    by = c("ORF_id_tr" = "orf_id")
  ) %>%
  
  # Remove irrelevant columns related to ORF detection criteria and features
  dplyr::select(-c(
    "start_dif", "stop_same", "orf_cds",
    "start_check", "transcript_biotype",
  )) %>%
  
  # Rename column `ORF_id_tr` to `orf_id` for consistency
  dplyr::rename(orf_id = ORF_id_tr) %>%
  
  # Mark ORFs as translated if detected in more than 'min_number_of_samples'
  dplyr::mutate(
    translated = ifelse(number_patient_samples > min_number_of_samples, TRUE, FALSE),
    # Assign "ATG" as the default start codon
    start_codon = "ATG", 
    # Add a source column for better tracking
    source = "ORFquant"
  )

#union(setdiff(colnames(orfquant_orfs), colnames(price_orfs)), setdiff(colnames(price_orfs), colnames(orfquant_orfs)))
####################################################################################
# Combine both ORF tables into a single dataset
orfs <- rbind(orfquant_orfs, price_orfs)



####################################################################################


## intORFs reannotation




####################################################################################

## CDS reannotation

#Check CDS regions that are very similar to called ORFs based on CDS loci and 
#annotated proteins.
#Check whether PRICE ORF calls overlap with ORFquant calls

orfs <- orfs %>%
  # Join the 'cds_check' table on 'orf_id' to bring in 'same_as_cds'
  #dplyr::left_join(cds_check[, c("orf_id", "same_as_cds")], by = "orf_id") %>%
  # Replace missing values and update 'orf_category_new' based on conditions
  dplyr::mutate(
    same_as_cds = tidyr::replace_na(same_as_cds, replace = 0),
    similarity_score = tidyr::replace_na(similarity_score, replace = 0),
    orf_category_new = ifelse(
      same_as_cds == TRUE, 
      "ORF_annotated",
      ifelse(similarity_score > min_similarity_score, "ORF_annotated", orf_category_new)
    )
  ) %>%
  # Remove the 'same_as_cds' column as it's no longer needed
  dplyr::select(-same_as_cds)


orfs <- orfs %>%
  # Change 'orf_category_new' to 'lncORF' if 'gene_biotype' is 'lncRNA'
  dplyr::mutate(
    orf_category_new = ifelse(gene_biotype == "lncRNA", "lncORF", orf_category_new)
  ) %>%
  # Rename 'source' column to 'orf_caller' for consistency
  dplyr::rename(orf_caller = source)


####################################################################################

#FINAL OUPUT -------------
#no further modifications are done on the output file from here
outfile="harmonised_orfs.csv"
write.table(orfs, file = outfile,
            sep = ",", 
            quote = F, 
            row.names = F)
