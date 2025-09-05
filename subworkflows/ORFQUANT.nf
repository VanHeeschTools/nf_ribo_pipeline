include { prepare_orfquant ; orfquant ; fix_orfquant } from '../modules/orfquant.nf'
include { write_collected_paths } from '../modules/helperFunctions.nf'

// Run ORFquant using merged riboseqc output file
workflow ORFQUANT {
    take:
    for_orfquant_files      // Path, RiboseQC output files
    orfquant_annotation     // Path, ORFquant Rannot file location
    orfquant_annot_package  // Path, BSgenome index directory
    package_install_loc     // Path, location where BSgenome package is installed
    pandoc_dir              // Path, location of directory where pandoc is located
    outdir                  // Path, output directory

    main:
    // Collect RiboseQC ouput paths into a single channel
    collected_paths = for_orfquant_files
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

    // Define ORFquant subworkflow output
    orfquant_orf_gtf = fix_orfquant.out.orfquant_gtf
    orfquant_orfs = orfquant.out.orfquant_orfs

    emit:
    orfquant_orf_gtf
    orfquant_orfs
}
