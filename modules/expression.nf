// Remove filtred out orf_ids from the orfcaller psites bed file
process filter_removed_orf_ids{

    label "filter_removed_orf_ids"
    publishDir "${outdir}/annotation", mode: 'copy'

    input:
    //TODO: make sure the removed_orf_ids are all removed from the orf_table
        path removed_orf_ids
        path orfcaller_psites
        val outdir

    output:
        path "combined_psites_filtered.bed", emit: orfcaller_psites_filtered

    script:
        """
        grep -v -w -F -f  ${removed_orf_ids} ${orfcaller_psites} > combined_psites_filtered.bed
        """
}

// Intersect the combined ORFcaller BED with p-site positions with p-sites from the samples
process intersect_psites {

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
        -s \
        -sorted > "${sample_id}_intersect.bed"
        """
}

// Create a matrix object for raw P-sites and PPM for each ORF and
// each sample included in the cohort
process ppm_matrix {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/orf_expression", mode: 'copy'

    input:
        path ref_psite_bed
        path sample_intersect_bed
        val outdir

    output:
        path "orf_table_psites_permillion.csv", emit: ppm_matrix
        path "orf_table_psites.csv", emit: psite_matrix

    script:
    """
        psite_matrix.R \
        "${ref_psite_bed}" \
        "${sample_intersect_bed}"
        """
}

// Add the expression information to the harmonised orf table
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
        # combined_ppm_table\$Total_number_samples <-  rowSums(combined_ppm_table[ , sample_cols] >= 1)

        combined_ppm_table\$Total_number_samples <- rowSums(
            combined_ppm_table[, sample_cols, drop = FALSE] >= 1
        )

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

// Create plot of translated Canonical and Non-canonical ORFs for MultiQC
process multiqc_expression_plot{

    label "Ribo_Seq_R_scripts"

    input:
        path harmonised_orf_table
        path ppm_matrix
        val outdir

    output:
        path "canonical_orf_counts_mqc.txt", emit: canonical_orf_counts

    script:
        """
        function_to_run="canonical_count"
        multiqc_tables.R \${function_to_run} "${harmonised_orf_table}" "${ppm_matrix}"
        """
}