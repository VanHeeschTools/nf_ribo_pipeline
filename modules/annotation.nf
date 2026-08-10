// Create reference rds file with correct cds start and stop 
process correct_reference_cds {

    label "Ribo_Seq_R"

    input:
        path orfcaller_gtf
        val type
        path reference_protein_fa
        path package_install_loc

    output:
        path "${orfcaller_gtf.baseName}_correct_cds.rds", emit: reference_cds_rds

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        process_ORF_and_Ref_gtf.R ${orfcaller_gtf} ${type} ${reference_protein_fa} ${package_install_loc}
        """
}

// Annotate ORFcaller output
process get_orf_category{

    label "Ribo_Seq_R"

    input:
        path orfcaller_gtf
        path reference_gtf         // Path, input reference gtf file
        path ref_cds_rds           // RDS file altered CDS regions that have a proper start and stop
        path package_install_loc   // Path, BSgenome package install location

    output:
        path "${orfcaller_gtf.baseName}_orfs.csv", emit: basic_orf_table

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        get_orf_categories.R \
        "${reference_gtf}" \
        "${orfcaller_gtf}" \
        "${ref_cds_rds}" \
        "${package_install_loc}" \
        "${orfcaller_gtf.baseName}" # ORFcaller name
        """
}

// Combine annotate_orfs results into a single coord sorted csv file
process harmonise_orfs {

    label "Ribo_Seq_R"
    input:
        val orfcaller_tables // List of paths to ORFcaller csv files
    
    output:
        path "harmonised_orf_table.csv",            emit: harmonised_orf_table
        path "removed_orf_ids.txt",                 emit: removed_orf_ids
        path "orf_protein_sequences.fa.gz",         emit: orf_protein_sequences
        path "orf_protein_sequences_M_start.fa.gz", emit: orf_protein_sequences_m_start
        path "orf_dna_sequences.fa.gz",             emit: orf_dna_sequences
        path "harmonised_orf_table.gtf",            emit: harmonised_orf_table_gtf
        path "orfcaller_orf_categories_mqc.txt",    emit: orfcaller_multiq
        path "merged_orf_categories_mqc.txt",       emit: merged_multiqc, optional: true
        path "merged_orf_caller_count_mqc.txt",     emit: caller_count_multiqc, optional: true

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        orf_harmonisation.R \
        ${orfcaller_tables.join(' ')}
        """
}

// Takes an existing ORF table df and converts it to a gtf-like file
process convert_table_to_gtf {

    label "Ribo_Seq_R"
    input:
        val orf_table
    
    output:
        path "harmonised_orf_table.gtf", emit: harmonised_orf_table_gtf

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        convert_to_gtf.R \
        ${orf_table}
        """
}