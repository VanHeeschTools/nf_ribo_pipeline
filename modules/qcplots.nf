// Create RiboseQC statistics plots for MultiQC
process riboseqc_plots {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/qc", mode: 'copy'

    input:
        path riboseqc_all
        val outdir

    output:
        path "Metagene_profile_combined_mqc.png", emit: metagene_plot
        path "Periodicity_bar_combined_mqc.png", emit: periodicity_plot

    script:
        """
        function_to_run="periodicity_plot"
        multiqc_tables.R \${function_to_run} ${riboseqc_all.join(' ')}
        """
}

// Create RiboseQC statistics tables for MultiQC
process riboseqc_tables {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/qc", mode: 'copy'

    input:
        val riboseqc_all
        val outdir

    output:
        path "inframe_percentages_mqc.txt", emit: riboseqc_inframe_percentages
        path "riboseqc_read_categories_counts_mqc.txt", emit: riboseqc_category_counts

    script:
        """
        function_to_run="riboseqc_tables"
        multiqc_tables.R \${function_to_run} ${riboseqc_all.join(' ')}
        """
}
