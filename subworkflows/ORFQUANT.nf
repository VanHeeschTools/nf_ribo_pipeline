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

    collected_paths = orfquant_psites
    .map { meta, path -> path }
    .collect()

    //write_collected_paths(collected_paths)

    prepare_orfquant(collected_paths,
                     orfquant_prefix,
                     outdir)

    orfquant(prepare_orfquant.out.psites_merged,
             orfquant_prefix,
             orfquant_annotation,
             pandoc_dir,
             orfquant_annot_package,
             package_install_loc,
             outdir)

    // Fixed the gtf output of orfquant which is mostly relevant if only orfquant is used
    fix_orfquant(orfquant.out.orfquant_orfs,
                 orfquant_annotation,
                 orfquant_prefix,
                 package_install_loc,
                 outdir)

    orfquant_orf_gtf = fix_orfquant.out.orfquant_gtf
    orfquant_orfs = orfquant.out.orfquant_orfs

    emit:
    orfquant_orf_gtf
    orfquant_orfs

}
