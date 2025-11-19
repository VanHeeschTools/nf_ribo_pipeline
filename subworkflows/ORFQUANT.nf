include { prepare_orfquant ; orfquant ; fix_orfquant } from '../modules/orfquant.nf'
// Run ORFquant using merged riboseqc output file
workflow ORFQUANT {
    take:
    for_orfquant_files      // Path, RiboseQC output files
    orfquant_annotation     // Path, ORFquant Rannot file location
    reference_gtf           // Path, input reference gtf file
    package_install_loc     // Path, location where BSgenome package is installed
    outdir                  // Path, output directory

    main:
    // Collect RiboseQC ouput paths into a single channel
    collected_paths = for_orfquant_files
        .map { _meta, path -> path }
        .collect()

    // Merge RiboseQC output 
    prepare_orfquant(
        collected_paths,
        outdir
    )

    // Run ORFquant using merged RiboseQC output
    orfquant(
        prepare_orfquant.out.psites_merged,
        orfquant_annotation,
        package_install_loc,
        outdir
    )

    // Corrects the IDs of the ORFquant GTF and adds plus three to the end coordinates of CDS
    fix_orfquant(
        orfquant.out.orfquant_orfs,
        orfquant_annotation,
        reference_gtf,
        package_install_loc,
        outdir
    )

    // Define ORFquant subworkflow output
    orfquant_orf_gtf = fix_orfquant.out.orfquant_gtf
    orfquant_orfs = orfquant.out.orfquant_orfs

    emit:
    orfquant_orf_gtf
    orfquant_orfs
}
