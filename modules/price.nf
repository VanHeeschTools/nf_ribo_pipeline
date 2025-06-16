process price_index {

    label "price"

    input:
    path fasta        // Genome fasta used for alingment
    path gtf          // Transcriptome GTF used for alignment
    val gedi_exec_loc // Location of gedi installation, until containerisation works

    output:
    path "PRICE_index.oml", emit: price_index
    path "${gtf.baseName}.*"
    path "${fasta.baseName}.*"

    script:
    """
    ${gedi_exec_loc}/gedi -e IndexGenome \
        -s "${fasta}" \
        -a "${gtf}" \
        -f "." \
        -o "PRICE_index.oml" \
        -nobowtie \
        -nostar \
        -nokallisto
    """

}

process merge_price_bams{
    label "samtools"
    publishDir "${outdir}/price", mode: 'copy'

    input:
    path bamlist
    path outdir

    output:
    path "star_end2end_merged_sorted.bam", emit: merged_end2end_bam

    script:
    """
    samtools merge -@ $task.cpus -b ${bamlist} star_end2end_merged.bam 
    samtools sort -@ $task.cpus -o star_end2end_merged_sorted.bam star_end2end_merged.bam
    samtools index star_end2end_merged_sorted.bam
    """

}

process price {

    label "price"
    publishDir "${outdir}/price", mode: 'copy'

    input:
    path bamlist      // File that lists all BAM files
    path price_index  // Index for PRICE
    val gedi_exec_loc // Location of gedi installation, until containerisation works
    val outdir        // Output directory

    output:
    tuple val("PRICE"), path("PRICE.orfs.cit.bed"), emit: price_orfs

    script:
    """
    ${gedi_exec_loc}/gedi -e Price \
        -reads ${bamlist} \
        -genomic ${price_index} \
        -prefix "PRICE"
    
    ${gedi_exec_loc}/gedi Nashorn -e \
        'load("'PRICE.orfs.cit'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' \
        > "PRICE.orfs.cit.bed"
    """

}

process price_to_gtf{
    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/price", mode: 'copy'

    input:
    tuple val(orfcaller), val(price_bed_file)
    val outdir

    output:
    path "PRICE.gtf", emit: price_gtf

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
