include { get_orf_category; annotate_orfs ; harmonise_orfs; merge_orfcaller_gtf } from "../modules/annotation.nf"

workflow ANNOTATION {
    take:
    reference_gtf           // Path, input gtf file
    package_install_loc     // Path, Location where BSgenome R package is installed
    run_ribotie             // Bool, True if RiboTIE should be run in the pipeline
    orf_gtf_bed             // Path, bed-like file of ORF and ref CDS overlap
    ref_cds_rds             // Path, RDS file of altered reference CDS
    outdir                  // Path, output directory

    main:
    // Load ORFcaller gtf and annotates the ORFs
    get_orf_category(
        orf_gtf_bed,
        reference_gtf,
        ref_cds_rds,
        package_install_loc,
        outdir
    )

    //TODO: put this logic dynamically in the orf_harmonisation_script
    if (run_ribotie) {
        // Collect ORF tables and sort them alphabetically
        harmonise_input = get_orf_category.out.basic_orf_table
            .collect()
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1], sorted[2])
            }
    }
    else {
        // Collect ORF tables and sort them alphabetically
        harmonise_input = get_orf_category.out.basic_orf_table
            .collect()
            .map { files ->
                def sorted = files.sort { it.name }
                tuple(sorted[0], sorted[1], null)
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
