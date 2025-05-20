process annotate_orfs {

    label "annotate_orfs"
    publishDir "${outdir}/annotate_orfs", mode: 'copy'

    input:
    tuple val(orfcaller), val (orfcaller_output) // Output of the orf caller (bed or ORFquant object)
    val orfcaller_psites
    val ref_psites
    val reference_gtf                    // gtf file
    val package_install_loc
    val orfquant_annot_package
    val outdir

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

process harmonise_orfs {

    label "annotate_orfs"
    publishDir "${outdir}/harmonise_orfs", mode: 'copy'

    input:
    tuple path(orfquant_table), path (price_table)
    val outdir

    output:
    path "harmonised_table.csv", emit: harmonised_orf_table
    path "removed_orf_ids.txt", emit: removed_orf_ids

    script:
    """
    orf_harmonisation.R \
    "${orfquant_table}" \
    "${price_table}" 
    """
}
