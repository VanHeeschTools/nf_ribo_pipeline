#!/usr/bin/env Rscript

# ============================================================
# FASTQ length counter for MultiQC lineplot (20–40 nt only)
# ============================================================
# Usage:
#   Rscript fastq_len.R input.fastq.gz output.txt
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

sample_id  <- args[1]
infile <- args[2]
outfile <- paste0(sample_id,"_size_distribution_mqc.txt")

count_lengths <- function(infile) {
  # Open input file (gzipped or plain text)
  con <- if (grepl("\\.gz$", infile)) {
    gzfile(infile, "rt")
  } else {
    file(infile, "rt")
  }
  
  # Preallocate a counter for read lengths 1–41
  # We will only use indices 20–40
  length_counts <- integer(41)
  
  repeat {
    # Read 40,000 lines = 10,000 FASTQ records (4 lines each)
    lines <- readLines(con, n = 40000)
    if (length(lines) == 0) break  # End of file
    
    # Extract sequence lines (2nd line in each 4-line FASTQ block)
    seq_lines <- lines[seq(2, length(lines), by = 4)]
    
    # Compute read lengths for this chunk
    read_lengths <- nchar(seq_lines)
    
    # Update counts only for lengths 20–40
    for (read_len in read_lengths) {
      if (read_len >= 20 && read_len <= 40) {
        length_counts[read_len] <- length_counts[read_len] + 1L
      }
    }
  }
  
  close(con)
  
  # Build result data.frame with only lengths that occurred
  result <- data.frame(
    length = 20:40,
    count  = length_counts[20:40]
  )
  
  result[result$count > 0, ]
}

# Run counting
len_df <- count_lengths(infile)

# Save output: two columns (length, count) without headers
write.table(
  len_df,
  file      = outfile,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)