include { orfcaller_psites ; reference_psites; sample_psites ; merge_orfcaller_psites;} from "../modules/process_bed.nf"
include { convert_table_to_gtf } from "../modules/annotation.nf"
// Obtain P0 sites of all reference transcripts and all predicted ORFs
// Required to classify intORFs and to calculate PPM
workflow PSITE {
    take:
    orfcaller_gtf        // Path, ORFcaller output in gtf format
    reference_gtf        // Path, input reference gtf
    reference_protein_fa // Path, input reference protein fasta
    package_install_loc  // Path, location where BSgenome package is installed
    run_quantify_existing
    existing_orf_table

    main:
    // If running quantification on existing ORF list, convert ORF list to gtf-like format
    if (run_quantify_existing){
        convert_table_to_gtf(
            existing_orf_table
        )
        orfcaller_gtf = convert_table_to_gtf.out
        ref_cds_rds = channel.empty() // Annotation won't be run in this mode so not required

    } else {
        // Create reference in-frame bed file for the reference gtf only when annotation is run
        reference_psites(
            reference_gtf,
            "transcript_id",
            reference_protein_fa,
            package_install_loc,
        )
        ref_cds_rds = reference_psites.out.reference_cds_rds
    }

    // Create reference in-frame bed file for the ORF caller
    orfcaller_psites(
        orfcaller_gtf,
        "ORF_id",
        reference_protein_fa,
        package_install_loc,
    )

    // Merge the bed files of all ORFcallers
    merge_orfcaller_psites(
        orfcaller_psites.out.orf_psite_bed.collect(),
    )

    // Define PSITE subworkflow output
    orfcaller_psites = merge_orfcaller_psites.out.combined_psites

    emit:
    // p0 locations in ORF callers
    orfcaller_psites
    // RDS file of the reference CDS
    ref_cds_rds
}
