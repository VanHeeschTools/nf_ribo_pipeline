include { prepare_orfquant; orfquant; fix_orfquant } from '../modules/orfquant.nf'
include { write_psites_paths } from '../modules/helperFunctions.nf'

workflow ORFQUANT {

    take:
    orfquant_psites
    orfquant_annotation
    orfquant_annot_package
    package_install_loc
    pandoc_dir
    orfquant_prefix
    outdir

    main:

    orfquant_psites
    .map { meta, path -> path } // Extract only the paths
    .collect()                   // Collect all paths into a single list
    .set { collected_paths }     // Set this list into a new channel
    write_psites_paths(collected_paths)

    prepare_orfquant(write_psites_paths.out.psites_file_channel,
                     orfquant_prefix,
                     outdir)

    orfquant(prepare_orfquant.out.psites_merged,
             orfquant_prefix,
             orfquant_annotation,
             pandoc_dir,
             orfquant_annot_package,
             package_install_loc,
             outdir
             )

    fix_orfquant(orfquant.out.orfquant_results_file,
                 orfquant_annotation,
                 orfquant_prefix)

    orfquant_orfs = fix_orfquant.out.orfquant_gtf
    orfquant_results_file = orfquant.out.orfquant_results_file

    emit:
    orfquant_orfs
    orfquant_results_file

}
