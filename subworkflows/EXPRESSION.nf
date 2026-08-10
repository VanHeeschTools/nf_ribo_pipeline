include { filter_removed_orf_ids; intersect_psites ; ppm_matrix ; multiqc_expression_plot; } from "../modules/expression.nf"
include { orfcaller_psites ; sample_psites ; merge_orfcaller_psites; } from "../modules/process_bed.nf"
include { convert_table_to_gtf } from "../modules/annotation.nf"

workflow EXPRESSION {
    take:
    orfcaller_gtf         // Path, ORFcaller output in gtf format
    package_install_loc   // Path, Location where BSgenome R package is installed
    for_orfquant_files    // Path, RiboseQC output files
    harmonised_orf_table  // Path, harmonised orf table csv file
    removed_orf_ids       // Path, txt file of filtered out ORF ids
    run_quantify_existing // Bool, True if quantification of samples needs to be done on existing ORF table
    existing_orf_table    // Path, existing ORF table from previous pipeline run

    main:

    // If running quantification on existing ORF list, convert ORF list to gtf-like format
    if (run_quantify_existing){
        convert_table_to_gtf(existing_orf_table)
        orfcaller_gtf = convert_table_to_gtf.out
    } 

    // Create bed file for the ORF caller showing all coordinates and frame information
    orfcaller_psites(
        orfcaller_gtf,
        "ORF_id",
        package_install_loc,
    )

    // Merge the bed files containing coords and frame information of all ORFcallers
    merge_orfcaller_psites(
        orfcaller_psites.out.orf_psite_bed.collect(),
    )

    // Remove filtered out ORF ids from the merged ORFcaller p0 site bed file
    // Will be skipped if running quantification on existing ORF list
    if (!run_quantify_existing){
        filter_removed_orf_ids(
            removed_orf_ids,
            merge_orfcaller_psites.out.combined_psites,
        )
        orfcaller_psites_filtered = filter_removed_orf_ids.out.orfcaller_psites_filtered
    } else {
        orfcaller_psites_filtered = merge_orfcaller_psites.out.combined_psites
    }

    // Create sample P-site files
    sample_psites(
        for_orfquant_files,
    )

    // Create intersect between P-sites and ORF locations to find all overlap
    intersect_psites(
        sample_psites.out.sample_psite_bed,
        orfcaller_psites_filtered,
    )

    // Calculate PPM matrices
    ppm_matrix(
        orfcaller_psites_filtered,
        intersect_psites.out.sample_intersect.collect(),
        harmonised_orf_table
    )

    // Define output ppm_matrix
    ppm_matrix_csv = ppm_matrix.out.ppm_matrix

    // Create MultiQC plot to show amount of expressed canonical and non-canonical ORFs
    multiqc_expression_plot(
        harmonised_orf_table,
        ppm_matrix_csv,
    )

    emit:
    // PPM matrix csv file
    ppm_matrix_csv
    // P-site matrix csv file
    psite_matrix_csv = ppm_matrix.out.psite_matrix
    // PIF and Uniformity score calculation
    orf_table_translation_scores = ppm_matrix.out.orf_table_translation_scores
    // EXPRESSION subworkflow statistics for MultiQC
    multiqc_expression_plot_txt = multiqc_expression_plot.out
}
