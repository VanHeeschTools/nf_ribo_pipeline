// Create RiboseQC statistics tables for MultiQC
process riboseqc_tables {

    label "Ribo_Seq_R_scripts"

    input:
        val riboseqc_all

    output:
        path "inframe_percentages_mqc.txt", emit: riboseqc_inframe_percentages
        path "riboseqc_read_categories_counts_mqc.txt", emit: riboseqc_category_counts

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        function_to_run="riboseqc_tables"
        multiqc_tables.R \${function_to_run} ${riboseqc_all.join(' ')}
        """
}
