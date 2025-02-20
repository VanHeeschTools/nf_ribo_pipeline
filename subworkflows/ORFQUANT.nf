include { prepare_orfquant; orfquant; fix_orfquant } from '../modules/orfquant.nf'
include { write_collected_paths } from '../modules/helperFunctions.nf'

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
    .map { meta, path -> path }
    .collect()
    .set { collected_paths }
    write_collected_paths(collected_paths)

    prepare_orfquant(write_collected_paths.out.collected_file_channel,
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
