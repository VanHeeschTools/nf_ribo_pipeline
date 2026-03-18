include { orfcaller_psites ; reference_psites; sample_psites ; merge_orfcaller_psites; orf_ref_p0_intersect } from "../modules/process_bed.nf"

// Obtain P0 sites of all reference transcripts and all predicted ORFs
// Required to classify intORFs and to calculate PPM
workflow PSITE {
    take:
    orfcaller_gtf        // Path, ORFcaller output in gtf format
    reference_gtf        // Path, input reference gtf
    reference_protein_fa // Path, input reference protein fasta
    package_install_loc  // Path, location where BSgenome package is installed

    main:
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

    // Create reference in-frame bed file for the reference gtf
    reference_psites(
        reference_gtf,
        "transcript_id",
        reference_protein_fa,
        package_install_loc,
    )

    // Define PSITE subworkflow output
    orfcaller_psites = merge_orfcaller_psites.out.combined_psites
    ref_cds_rds = reference_psites.out.reference_cds_rds

    emit:
    orfcaller_psites
    ref_cds_rds
}
