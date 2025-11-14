#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

riboseqc_file <- args[1]

# Set location of installed packages for this R version

suppressPackageStartupMessages({
    library(RiboseQC)
    #library(data.table)
})

# Get basename of the file
fname <- gsub(pattern = "_for_ORFquant", x = basename(riboseqc_file), replacement = "")
print(fname)

# Load the file and capture the name of the loaded object
loaded_obj_names <- load(riboseqc_file)
print(loaded_obj_names)

# Assume only one object is loaded; get that object
my_obj <- get(loaded_obj_names[1])

# Extract p-sites from list of objects
p_sites <- data.frame(my_obj$P_sites_uniq)

if(nrow(p_sites) == 0) {
    stop("No p-sites found in the loaded object.")
}

# Extract columns to create BED file
bed <- data.frame(
    chrom = p_sites$seqnames,
    chromStart = as.integer(p_sites$start - 1), # Minus 1 is used to create proper bed file 
    chromEnd = as.integer(p_sites$end),
    name = ".",
    score = p_sites$score,
    strand = p_sites$strand
)


write.table(bed, file = paste0(fname, "_psites.bed"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
