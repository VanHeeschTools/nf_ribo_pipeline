process riboseqc_plots {

    label "riboseqc_plots"
    publishDir "${outdir}/qc", mode: 'copy'

    input:
    val data_files
    val outdir
    val pandoc_dir
    val render_file
    val orfquant_prefix

    output:
    "*.html"

    script:
    """
    Rscript riboseqc_html.R \
    ${input_files}
    ${pandoc_dir} \
    ${render_file} \
    ${orfquant_prefix}
    """
}

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
