include { annotate_orfs ; harmonise_orfs; convert_csv_to_gtf } from "../modules/annotation.nf"

workflow ANNOTATION {
    take:
    orfcaller_output        // Path, ORFcaller output file
    orfcaller_psites        // Path, merged p0 sites of all used ORFcallers
    ref_psites              // Path, p0 sites of reference gtf
    reference_gtf           // Path, input gtf file
    package_install_loc     // Path, Location where BSgenome R package is installed
    orfquant_annot_package  // Path, BSgenome index directory
    run_ribotie             // Bool, True if RiboTIE should be run in the pipeline
    orfcaller_gtf
    outdir                  // Path, output directory

    main:
    // Parses and annotates the output of the ORFcallers
    annotate_orfs(
        orfcaller_output,
        orfcaller_psites,
        ref_psites,
        reference_gtf,
        package_install_loc,
        orfquant_annot_package,
        outdir
    )

    // Collect ORF tables and sort them alphabetically
    if (run_ribotie) {
        harmonise_input = annotate_orfs.out.basic_orf_table
            .collect()
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1], sorted[2])
            }
        orfcaller_gtf_sorted = orfcaller_gtf
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1], sorted[2])
            }
    }
    else {
        harmonise_input = annotate_orfs.out.basic_orf_table
            .collect()
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1], null)
            }
        orfcaller_gtf_sorted = orfcaller_gtf
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1])
            }
    }

    // Combines the ORFcaller annotated csv files
    harmonise_orfs(
        harmonise_input,
        outdir,
    )

    // Define subworkflow output
    harmonised_orf_table = harmonise_orfs.out.harmonised_orf_table
    removed_orf_ids = harmonise_orfs.out.removed_orf_ids

    convert_csv_to_gtf(harmonised_orf_table, 
                    orfcaller_gtf_sorted,
                    outdir)


    // Annotation multiqc output files
    orfcaller_multiq = harmonise_orfs.out.orfcaller_multiq
    merged_multiqc = harmonise_orfs.out.merged_multiqc
    caller_count_multiqc = harmonise_orfs.out.caller_count_multiqc
    annotation_multiqc = orfcaller_multiq.mix(merged_multiqc, caller_count_multiqc)

    emit:
    harmonised_orf_table
    removed_orf_ids
    annotation_multiqc
}
