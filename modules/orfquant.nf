// Merge RiboseQC output 
process prepare_orfquant {

    label "Ribo_Seq_R"

    input:
        val collected_paths

    output:
        path "Merged_for_ORFquant", emit: psites_merged
        path "file_paths.txt",      emit: for_orfquant_file_paths_txt

    when:
        task.ext.when == null || task.ext.when

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

    label "Ribo_Seq_R"

    input:
        val psites_merged
        val rannot
        val package_install_loc

    output:
        path "output_final_ORFquant_results", emit: orfquant_orfs

    when:
        task.ext.when == null || task.ext.when

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
    label "Ribo_Seq_R"

    input:
        path(orfquant_orfs)
        path rannot
        path reference_gtf
        val package_install_loc

    output:
        path "ORFquant.gtf", emit: orfquant_gtf

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        fix_orfquant_output.R \
            ${orfquant_orfs} \
            ${rannot} \
            ${package_install_loc} \
            ${reference_gtf}
        """
}
