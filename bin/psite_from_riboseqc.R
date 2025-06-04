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

# ORFquant files are R object files
#load(riboseqc_file)

# Load the file and capture the name of the loaded object
loaded_obj_names <- load(riboseqc_file)
print(loaded_obj_names)

# Assume only one object is loaded; get that object
my_obj <- get(loaded_obj_names[1])

# Extract p-sites from list of objects
# TODO: test with P_sites_all
p_sites <- data.frame(my_obj$P_sites_uniq)

if(nrow(p_sites) == 0) {
    stop("No p-sites found in the loaded object.")
}

# Extract columns to create BED file
bed <- data.frame(p_sites$seqnames, p_sites$start, p_sites$end, ".", p_sites$score, p_sites$strand)
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

#data.table::fwrite(bed, paste0(fname, "_psites.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
write.table(bed, file = paste0(fname, "_psites.bed"), quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
