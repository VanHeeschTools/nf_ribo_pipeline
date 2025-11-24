// Annotate ORFcaller output
process get_orf_category{

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/annotate_orfs", mode: 'copy'

    input:
    tuple path(orfcaller_gtf), path(cds_orf_bed_file)
    path reference_gtf         // Path, input reference gtf file
    path ref_cds_rds           // RDS file altered CDS regions that have a proper start and stop
    path package_install_loc   // Path, BSgenome package install location
    val outdir                 // Path, output directory

    output:
    path "${orfcaller_gtf.baseName}_orfs.csv", emit: basic_orf_table

    script:
    """
    get_orf_categories.R \
    "${reference_gtf}" \
    "${orfcaller_gtf}" \
    "${cds_orf_bed_file}" \
    "${ref_cds_rds}" \
    "${package_install_loc}" \
    "${orfcaller_gtf.baseName}" # ORFcaller name
    """
}

// Combine annotate_orfs results into a single coord sorted csv file
process harmonise_orfs {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/harmonise_orfs", mode: 'copy'
    publishDir "${outdir}/final_orf_table", mode: 'copy', pattern: 'orf_sequences.fa.gz'
    publishDir "${outdir}/final_orf_table", mode: 'copy', pattern: 'orf_harmonised.gtf'

    input:
    val orfcaller_tables
    val outdir                                                        
    
    output:
    path "harmonised_orf_table.csv", emit: harmonised_orf_table
    path "removed_orf_ids.txt", emit: removed_orf_ids
    path "orf_protein_sequences.fa.gz"
    path "orf_dna_sequences.fa.gz"
    path "harmonised_orf_table.gtf"

    path "orfcaller_orf_categories_mqc.txt", emit: orfcaller_multiq
    path "merged_orf_categories_mqc.txt", emit: merged_multiqc
    path "merged_orf_caller_count_mqc.txt", emit: caller_count_multiqc

    script:
    """
    orf_harmonisation.R \
    ${orfcaller_tables.join(' ')}
    """
}

