include { ref_psites; sample_psites} from "../modules/process_bed.nf"
include { intersect_psites; ppm_matrix } from "../modules/expression.nf"

workflow EXPRESSION {

    take:
    orfs
    riboseqc_results
    package_install_loc
    outdir

    main:

    // Create sample P-site files
    sample_psites(riboseqc_results,
                  package_install_loc,
                  outdir
    )
    
    // Create reference in-frame bed file for the ORF caller
    ref_psites(orfs,
               outdir
    )


    intersect_input = sample_psites.out.sample_psite_bed.cross(ref_psites.out.ref_psite_bed)

    // Create intersect between P-sites and in-frame ORF locations
    // TODO: should run for each sample, but only runs for the first sample now
    intersect_psites(sample_psites.out.sample_psite_bed,
                     ref_psites.out.ref_psite_bed.first(),
                     outdir
    )

    collected_paths = intersect_psites.out.sample_intersect_bed
    .map { meta, path -> path }
    .collect()
    //write_collected_paths(collected_paths)
    //TODO: check if this format is correct for ppm_matrix

    // Calculate PPM matrices
    ppm_matrix(ref_psites.out.ref_psite_bed.first(),
               collected_paths,
               "orfcaller_name",
               outdir
    )

    ref_psite_bed = ref_psites.out.ref_psite_bed
    ppm_matrix_out = ppm_matrix.out.ppm_matrix
    raw_matrix_out = ppm_matrix.out.psite_matrix

    emit:
    ref_psite_bed
    ppm_matrix_out
    raw_matrix_out

}
