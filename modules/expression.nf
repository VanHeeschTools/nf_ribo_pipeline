process intersect_psites {

    // Intersect a reference BED with p-site positions with p-sites from a sample

    tag "${sample_id}"
    label "intersect_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_psite_bed)
    path ref_psite_bed
    val outdir

    output:
    path "${sample_id}_intersect.bed", emit: sample_intersect

    script:
    """
    bedtools intersect \
      -a ${sample_psite_bed} \
      -b ${ref_psite_bed} \
      -wa \
      -wb \
      -header \
      -f 1.00 \
      -s \
      -sorted > "${sample_id}_intersect.bed"
      """
}


process ppm_matrix {

    // Create a matrix object for raw P-sites and PPM for each ORF and
    // each sample included in the cohort

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/orf_expression", mode: 'copy'

    input:
    path ref_psite_bed
    path sample_intersect_bed

    val outdir

    output:
    path("orf_table_psites_permillion.csv"), emit: ppm_matrix
    path("orf_table_psites.csv"), emit: psite_matrix

    script:
    """
    psite_matrix.R \
    "${ref_psite_bed}" \
    "${sample_intersect_bed}"
    """
}

process expression_table{
    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/final_orf_table", mode: 'copy'


    input:
    val harmonised_orf_table
    val ppm_matrix
    val outdir

    output:
    path("final_orf_table.csv"), emit: final_orf_table

    script:
    """
    #!/usr/bin/env Rscript

    library(dplyr)

    # Load ORF table and PPM table
    orf_table <- read.delim("${harmonised_orf_table}", sep = ",") 
    combined_ppm_table <- read.csv("${ppm_matrix}", header = TRUE, row.names = NULL)

    # Calculate amount of samples with a PPM of 1 or higher
    sample_cols <- setdiff(names(combined_ppm_table), "orf_id")
    combined_ppm_table\$Total_number_samples <-  rowSums(combined_ppm_table[ , sample_cols] >= 1)

    # Join ORF table with PPM results
    orf_table_joined <- orf_table %>%
    left_join(
        combined_ppm_table %>%
        select(orf_id, Total_number_samples),
        by = "orf_id"
    )

    write.table(orf_table_joined, file = "final_orf_table.csv",
                sep = ",",
                quote = F,
                row.names = F)
    """
}