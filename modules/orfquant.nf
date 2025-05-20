process prepare_orfquant {

    label "orfquant_prep"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val collected_paths
    val orfquant_prefix
    val outdir

    output:
    path "${orfquant_prefix}_for_ORFquant", emit: psites_merged

    script:

    """
    printf "%s\n" "${collected_paths.join('\n')}" > file_paths.txt

    merge_psites.R \
    "file_paths.txt" \
    ${orfquant_prefix}
    """
}

process orfquant {

    label "orfquant"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val psites_merged
    val orfquant_prefix
    val rannot
    val pandoc_dir
    val orfquant_annot_package
    val package_install_loc
    val outdir

    output:
    tuple val("ORFquant"), path("${orfquant_prefix}_final_ORFquant_results"), emit: orfquant_orfs
    //path "${orfquant_prefix}_*"

    script:
    """
    run_ORFquant.R \
    ${psites_merged} \
    ${orfquant_prefix} \
    ${rannot} \
    $task.cpus \
    ${pandoc_dir} \
    ${orfquant_annot_package} \
    ${package_install_loc} 
    """
}

process fix_orfquant {

    // Fixes ORFquant GTF which has incorrect names and doesn't include the stop codon in the coords

    label "fix_orfquant"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    tuple val(orfcaller), path(orfquant_orfs)
    path rannot
    val orfquant_prefix
    val package_install_loc
    val outdir

    output:
    path "${orfquant_prefix}_Detected_ORFs_fixed.gtf", emit: orfquant_gtf
    //path "${orfquant_prefix}_Protein_sequences_fixed.fasta", emit: orfquant_fasta

    script:
    """
    fix_orfquant_output.R \
    ${orfquant_orfs} \
    ${rannot} \
    ${orfquant_prefix} \
    ${package_install_loc}
    """
}
