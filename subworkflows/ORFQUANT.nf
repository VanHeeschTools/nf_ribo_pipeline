include { prepare_orfquant ; orfquant ; fix_orfquant } from '../modules/orfquant.nf'
include { write_collected_paths } from '../modules/helperFunctions.nf'

// Run ORFquant using merged riboseqc output file
workflow ORFQUANT {
    take:
    orfquant_psites
    orfquant_annotation
    orfquant_annot_package
    package_install_loc
    pandoc_dir
    outdir

    main:

    // Collect RiboseQC ouput paths into a single channel
    collected_paths = orfquant_psites
        .map { _meta, path -> path }
        .collect()

    // Merge RiboseQC output 
    prepare_orfquant(
        collected_paths,
        outdir,
    )

    // Run ORFquant using merged RiboseQC output
    orfquant(
        prepare_orfquant.out.psites_merged,
        orfquant_annotation,
        pandoc_dir,
        orfquant_annot_package,
        package_install_loc,
        outdir,
    )

    // Corrects the IDs of the ORFquant GTF and adds plus three to the end coordinates of CDS
    fix_orfquant(
        orfquant.out.orfquant_orfs,
        orfquant_annotation,
        package_install_loc,
        outdir,
    )

    orfquant_orf_gtf = fix_orfquant.out.orfquant_gtf
    orfquant_orfs = orfquant.out.orfquant_orfs

    emit:
    orfquant_orf_gtf
    orfquant_orfs
}
