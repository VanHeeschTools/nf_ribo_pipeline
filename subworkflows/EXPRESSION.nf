include { sample_psites } from "../modules/process_bed.nf"
include { filter_removed_orf_ids ; intersect_psites ; ppm_matrix ; expression_table; multiqc_expression_plot } from "../modules/expression.nf"

workflow EXPRESSION {
    take:
    for_orfquant_files   // Path, RiboseQC output files
    harmonised_orf_table // Path, harmonised orf table csv file
    removed_orf_ids      // Path, txt file of filtered out ORF ids
    orfcaller_psites     // Path, merged p0 sites of all used ORFcallers
    outdir               // Path, output directory

    main:
    // Remove filtered out ORF ids from the merged ORFcaller p0 site bed file
    filter_removed_orf_ids(
        removed_orf_ids,
        orfcaller_psites
    )
    orfcaller_psites_filtered = filter_removed_orf_ids.out.orfcaller_psites_filtered

    // Create sample P-site files
    sample_psites(
        for_orfquant_files,
        outdir
    )

    // Create intersect between P-sites and ORF locations
    intersect_psites(
        sample_psites.out.sample_psite_bed,
        orfcaller_psites_filtered,
        outdir
    )

    // Obtain all intersect files before continuing with next step
    intersect_paths = intersect_psites.out.sample_intersect.collect()

    // Calculate PPM matrices
    ppm_matrix(
        orfcaller_psites_filtered,
        intersect_paths,
        outdir
    )

    // Define output ppm_matrix
    ppm_matrix = ppm_matrix.out.ppm_matrix

    // Merge the PPM results with the ORF table, creating the final output table
    expression_table(
        harmonised_orf_table,
        ppm_matrix,
        outdir
    )

    multiqc_expression_plot(
        harmonised_orf_table,
        ppm_matrix,
        outdir
    )

    multiqc_expression_plot_txt = multiqc_expression_plot.out


    emit:
    ppm_matrix
    multiqc_expression_plot_txt
}
