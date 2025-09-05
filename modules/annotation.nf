// Annotate ORFcaller output
process annotate_orfs {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/annotate_orfs", mode: 'copy'

    input:
    tuple val(orfcaller), val (orfcaller_output) // Output of the orf caller (bed or ORFquant object)
    val orfcaller_psites                         // P0 sites of orfcaller
    val ref_psites                               // P0 sites of reference transcripts
    val reference_gtf                            // Input gtf file
    val package_install_loc                      // Package install location
    val orfquant_annot_package                   // BSgenome location
    val outdir                                   // Path to output directory

    output:
    path "${orfcaller}_orfs.csv" , emit: basic_orf_table

    script:
    """
    orf_annotate.R \
    "${orfcaller_output}" \
    "${orfcaller_psites}" \
    "${ref_psites}" \
    "${reference_gtf}" \
    "${orfcaller}" \
    "${package_install_loc}" \
    "${orfquant_annot_package}" 
    """
}

// Combine annotate_orfs results into a single coord sorted csv file
process harmonise_orfs {

    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/harmonise_orfs", mode: 'copy'

    input:
    tuple path(orfquant_table), path(price_table), val(ribotie_table) // Tuple, path of orfquant and price annotation tables
    val outdir                                                         // Path to output directory

    output:
    path "harmonised_table.csv", emit: harmonised_orf_table
    path "removed_orf_ids.txt", emit: removed_orf_ids
    path "unfiltered_harmonised_table.csv"
    path "orfcaller_orf_categories_mqc.txt", emit: orfcaller_multiq
    path "merged_orf_categories_mqc.txt", emit: merged_multiqc
    path "merged_orf_caller_count_mqc.txt", emit: caller_count_multiqc

    //TODO: Add a proper preference checker in case of duplicate ORFs between ORFcallers
    script:
    """
    orf_harmonisation.R \
    "${orfquant_table}" \
    "${price_table}" \
    "${ribotie_table}"
    """
}
