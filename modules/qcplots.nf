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
    periodicity_plot.R ${riboseqc_all.join(' ')}
    """
}

// Create RiboseQC statistics tables for MultiQC
process riboseqc_tables {

    label "Ribo_Seq_R_scripts"

    input:
    val riboseqc_all

    output:
    path "riboseqc_frame_29nt_mqc.txt", emit: riboseqc_inframe_29
    path "riboseqc_read_categories_counts_mqc.txt", emit: riboseqc_category_counts

    script:

    """
    riboseqc_tables.R ${riboseqc_all.join(' ')}
    """
}
