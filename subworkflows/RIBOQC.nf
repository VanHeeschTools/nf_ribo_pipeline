include { riboseqc        } from '../modules/riboseqc.nf'
include { riboseqc_tables } from '../modules/qcplots.nf'

workflow RIBOQC {

    take:
    orfquant_annotation        // Path, ORFquant annotation file
    orfquant_annot_package     // Path, ORFquant annotation R package
    package_install_loc        // Path, location where BSgenome package is installed
    pandoc_dir                 // Path, location of pandoc for R HTML creation
    orfquant_bams              // List, output from ALIGNMENT subworkflow
    outdir                     // Path, output directory

    main:
    // Create riboseqc files
    riboseqc(orfquant_bams,
            outdir,
            orfquant_annotation,
            pandoc_dir,
            orfquant_annot_package,
            package_install_loc)

    // Create riboseqc tables for MultiQC
    riboseqc_tables(riboseqc.out.riboseqc_all.collect())
    riboseqc_inframe_29 = riboseqc_tables.out.riboseqc_inframe_29
    riboseqc_category_counts = riboseqc_tables.out.riboseqc_category_counts

    // Combine into one channel for MultiQC
    multiqc_riboseq = riboseqc_inframe_29.mix(riboseqc_category_counts)

    // Obtain ORFquant input files
    // Collect is done in the RIBOSEQ workflow
    for_orfquant_files = riboseqc.out.orfquant_psites

    emit:
    orfquant_annotation    // Used R annotation
    orfquant_annot_package // Used R package
    for_orfquant_files     // Files for ORFquant
    multiqc_riboseq        // Files for MultiQC

}
