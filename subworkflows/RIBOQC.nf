include { riboseqc; create_annotation } from "../modules/riboseqc.nf"
// include { multiqc } from "../modules/multiqc.nf"
include { riboseqc_tables } from '../modules/qcplots.nf'

workflow RIBOQC {

    take:
    orfquant_annotation        // ORFquant annotation file
    orfquant_annot_package     // ORFquant annotation R package
    package_install_loc        // Location where R package is installed
    pandoc_dir                 // Location of pandoc for R HTML creation
    orfquant_bams              // Output from ALIGNMENT subworkflow
    orfquant_annotation_exists // Boolean: whether to generate annotation
    gtf                        // Transcriptome used for STAR
    twobit                     // UCSC file format for the fasta
    outdir                     // Output directory

    //contaminants // RPF filtering statistics
    //star_output // STAR statistics per file
    //samtools_output // SAMTOOLS statistics per file

    main:
    //TODO: Fix this
    if(!orfquant_annotation_exists) {
        // Create riboseqc annotation
        create_annotation(gtf,
                          twobit,
                          package_install_loc,
                          orfquant_prefix)
        rannot_ch = create_annotation.out.orfquant_annotation
        package_ch = create_annotation.out.annotation_package
    } else {
        rannot_ch = orfquant_annotation
        package_ch = orfquant_annot_package
    }

    // Create riboseqc files
    riboseqc(orfquant_bams,
             outdir,
             rannot_ch,
             pandoc_dir,
             package_ch,
             package_install_loc)

    // Create riboseqc tables for MultiQC
    riboseqc_tables(riboseqc.out.riboseqc_all.collect())
  
    riboseqc_inframe_29 = riboseqc_tables.out.riboseqc_inframe_29
    riboseqc_category_counts = riboseqc_tables.out.riboseqc_category_counts
 
    multiqc_riboseq = riboseqc_inframe_29.mix(riboseqc_category_counts)

    for_orfquant_files = riboseqc.out.orfquant_psites

    emit:
    rannot_ch          // Used R annotation
    package_ch         // Used R package
    for_orfquant_files // Files for ORFquant
    multiqc_riboseq

}
