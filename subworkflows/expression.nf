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

    // Create intersect between P-sites and in-frame ORF locations
    intersect_psites(ref_psites.out.ref_bed,
                     sample_psites.out.samples_bed,
                     outdir
    )

    intersect_psites.out.intersect_bed
    .map { meta, path -> path }
    .collect()
    .set { collected_paths }
    //write_collected_paths(collected_paths)

    // Calculate PPM matrices
    ppm_matrix(ref_psites.out.ref_bed,
               "orfcaller_name",
               collected_paths,
               outdir
    )

    ref_bed = ref_psites.out.ref_bed
    ppm_matrix_out = ppm_matrix.out.ppm_matrix
    raw_matrix_out = ppm_matrix.out.psite_matrix

    emit:
    ref_bed
    ppm_matrix_out
    raw_matrix_out

}
