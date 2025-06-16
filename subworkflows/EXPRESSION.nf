include { sample_psites } from "../modules/process_bed.nf"
include { intersect_psites; ppm_matrix; expression_table } from "../modules/expression.nf"

workflow EXPRESSION {

    take:
    for_orfquant_files
    harmonised_orf_table
    removed_orf_ids
    orfcaller_psites
    outdir

    main:

    // Create sample P-site files
    sample_psites(for_orfquant_files,
                  outdir
    )
    
    // Create intersect between P-sites and in-frame ORF locations
    intersect_psites(sample_psites.out.sample_psite_bed,
                     orfcaller_psites,  //.first(),
                     outdir
    )

    intersect_paths = intersect_psites.out.sample_intersect.collect()
 
    // Calculate PPM matrices
    ppm_matrix(orfcaller_psites,
               intersect_paths,
               outdir
    )

    ppm_matrix = ppm_matrix.out.ppm_matrix

    // Merge the PPM results with the ORF table
    expression_table(harmonised_orf_table,
                     ppm_matrix,
                     outdir
    )

    emit:
    ppm_matrix


}
