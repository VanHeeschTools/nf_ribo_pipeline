include { sample_psites } from "../modules/process_bed.nf"
include { filter_removed_orf_ids ; intersect_psites ; ppm_matrix ; expression_table } from "../modules/expression.nf"

workflow EXPRESSION {
    take:
    for_orfquant_files
    harmonised_orf_table
    removed_orf_ids
    orfcaller_psites
    outdir

    main:
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

    intersect_paths = intersect_psites.out.sample_intersect.collect()

    // Calculate PPM matrices
    ppm_matrix(
        orfcaller_psites_filtered,
        intersect_paths,
        outdir
    )

    ppm_matrix = ppm_matrix.out.ppm_matrix

    // Merge the PPM results with the ORF table, creating the final output table
    expression_table(
        harmonised_orf_table,
        ppm_matrix,
        outdir
    )

    emit:
    ppm_matrix
}
