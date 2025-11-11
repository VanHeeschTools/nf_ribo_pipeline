// Merge RiboseQC output 
process prepare_orfquant {

    label "orfquant"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val collected_paths
    val outdir

    output:
    path "Merged_for_ORFquant", emit: psites_merged
    path "file_paths.txt"

    script:

    """
    #Collect all RiboseQC output paths in one file
    printf "%s\n" "${collected_paths.join('\n')}" > file_paths.txt

    merge_psites.R \
        "file_paths.txt"
    """
}

// Run ORFquant on merged psites level
process orfquant {

    label "orfquant"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val psites_merged
    val rannot
    val package_install_loc
    val outdir

    output:
    path "output_final_ORFquant_results", emit: orfquant_orfs

    script:
    """
    run_ORFquant.R \
        ${psites_merged} \
        ${rannot} \
        $task.cpus \
        ${package_install_loc} 
    """
}

// Fixes ORFquant GTF which has incorrect names and doesn't include the stop codon in the coords
process fix_orfquant {
    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    path(orfquant_orfs)
    path rannot
    path reference_gtf
    val package_install_loc
    val outdir

    output:
    path "ORFquant.gtf", emit: orfquant_gtf

    script:
    """
    fix_orfquant_output.R \
        ${orfquant_orfs} \
        ${rannot} \
        ${package_install_loc} \
        ${reference_gtf}
    """
}
