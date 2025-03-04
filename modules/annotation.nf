process annotate_orfs {

    label "annotate_orfs"
    publishDir "${outdir}/annotate_orfs", mode: 'copy'

    input:
    val orfs_loc                  // Output of the orf caller (bed or ORFquant object)
    val reference_gtf             // gtf file
    val orfcaller                 // Val, name of the orfcaller used
    val annotation_provider       // Val, name of the location where the IDs come from, ensembl or gencode
    val gencode_uniprot_file      // File, uniprot file required when annotation_provider is gencode
    val uniprot_protein_fasta_loc //File, location of uniprot fasta file
    val package_install_loc
    val orfquant_annot_package

    output:
    path "${orf_caller}_orfs.csv"

    script:
    """
    orf_annotate.R \
    "${orfs_loc}" \
    "${reference_gtf}" \
    "${orfcaller}" \
    "${annotation_provider}" \
    "${gencode_uniprot_file}" \
    "${uniprot_protein_fasta_loc}" \
    "${package_install_loc}" \
    "${orfquant_annot_package}" \
    $task.cpus
    """
}

process harmonise_orfs{

    // Should only be run when both ORFquant and PRICE are run

    label "harmonise_ofs"
    publishDir "${outdir}/harmonise_orfs", mode: 'copy'

    input:

    output:

    script:
    """
    orf_harmonisation.R \
    ""\
    ""
    # Needs a gtf, txdb
    # Needs the annotation r script output for price and orfquant

    """

}