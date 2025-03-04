include { ref_psites; sample_psites} from "../modules/process_bed.nf"
include { intersect_psites; ppm_matrix } from "../modules/expression.nf"

workflow EXPRESSION {

    take:
    orfs
    riboseqc_results
    package_install_loc
    orfcaller
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
    intersect_psites(sample_psites.out.sample_psite_bed,
                     ref_psites.out.ref_psite_bed.first(),
                     orfcaller,
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
               orfcaller,
               outdir
    )

    ref_psite_bed = ref_psites.out.ref_psite_bed
    ppm_matrix = ppm_matrix.out.ppm_matrix

    emit:
    ref_psite_bed
    ppm_matrix


}
