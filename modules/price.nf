// Create PRICE index file
process price_index {

    label "Ribo_Seq_tools"

    input:
        path fasta        // Genome fasta used for alingment
        path gtf          // Transcriptome GTF used for alignment

    output:
        path "PRICE_index.oml", emit: price_index

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        gedi -e IndexGenome \
            -s "${fasta}" \
            -a "${gtf}" \
            -f "." \
            -o "PRICE_index.oml" \
            -nobowtie \
            -nostar \
            -nokallisto
        """
}

// Merge all end2end BAM files into a single BAM to be used by PRICE
process merge_price_bams{
    label "Ribo_Seq_tools"

    input:
        path bam_files // File that lists all BAM files

    output:
        tuple path("star_end2end_merged_sorted.bam"), path("star_end2end_merged_sorted.bam.bai"), emit: merged_end2end_bam

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        samtools merge -@ $task.cpus star_end2end_merged.bam ${bam_files.join(' ')}
        samtools sort -@ $task.cpus -o star_end2end_merged_sorted.bam star_end2end_merged.bam
        samtools index star_end2end_merged_sorted.bam
        """
}

// Run PRICE on merged bam
process price {

    label "Ribo_Seq_tools"

    input:
        tuple path(merged_bam), path(merged_bam_bai)
        path price_index  // Index for PRICE

    output:
        path "PRICE.orfs.cit.bed", emit: price_orfs

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        export _JAVA_OPTIONS="${task.ext.java_opts}"
        gedi -e Price \
            -reads ${merged_bam} \
            -genomic ${price_index} \
            -prefix "PRICE"
        
        gedi Nashorn -e \
            'load("'PRICE.orfs.cit'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' \
            > "PRICE.orfs.cit.bed"
        """
}

// Convert PRICE output bed file to a semi gtf format only keeping the CDS rows
process price_to_gtf{
    label "Ribo_Seq_R"

    input:
        val price_bed_file

    output:
        path "PRICE.gtf", emit: price_gtf

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        #!/usr/bin/env Rscript

        library(rtracklayer)
        library(dplyr)
        library(stringr)
        # Load the BED file and generate ORF blocks
        price_orfs <- import.bed("${price_bed_file}")
        orf_ranges_list <- blocks(price_orfs)

        # Flatten the blocks and retain the original 'name' for each block
        orf_ranges <- unlist(orf_ranges_list, use.names = FALSE)
        orf_ranges\$ORF_id <- rep(price_orfs\$name, elementNROWS(orf_ranges_list))  # Keep 'name' as 'ORF_id'

        # Add additional columns
        orf_ranges\$type <- "CDS"
        orf_ranges\$source <- "PRICE"
        orf_ranges\$score <- "."
        orf_ranges\$frame <- "."

        # Extract gene and transcript info from the 'name' column
        orf_ranges_df <- as.data.frame(orf_ranges) %>%
        mutate(
            gene_id = str_split_i(ORF_id, "__", 1),
            transcript_id = str_split_i(ORF_id, "__", 2),
            start_codon = str_split_i(ORF_id, "__", 5)
        ) %>%
        dplyr::select(seqnames, start, end, score, strand, type, source, gene_id, transcript_id, ORF_id, start_codon)
        export.gff(orf_ranges_df, con = "PRICE.gtf")
        """
}
