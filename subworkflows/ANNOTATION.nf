include { get_orf_category; harmonise_orfs; convert_table_to_gtf; correct_reference_cds } from "../modules/annotation.nf"

workflow ANNOTATION {
    take:
    reference_gtf           // Path, input gtf file
    package_install_loc     // Path, Location where BSgenome R package is installed
    orfcaller_gtf           // Path, ORFcaller output in gtf format
    reference_protein_fa    // Path, input reference protein fasta

    main:
    // Create reference rds file with corrected cds start and stop location
    correct_reference_cds(
        reference_gtf,
        "transcript_id",
        reference_protein_fa,
        package_install_loc,
    )

    // Load ORFcaller gtf and annotates the ORFs
    get_orf_category(
        orfcaller_gtf,
        reference_gtf,
        correct_reference_cds.out.reference_cds_rds,
        package_install_loc,
    )

    // Collect the annotated ORF tables of all ORFcallers for harmonisation
    annotated_orf_tables = get_orf_category.out.basic_orf_table.collect()

    // Combines the ORFcaller annotated csv files into one harmonised ORF table
    harmonise_orfs(
        annotated_orf_tables
    )

    // Annotation multiqc output files
    orfcaller_multiq = harmonise_orfs.out.orfcaller_multiq
    merged_multiqc = harmonise_orfs.out.merged_multiqc
    caller_count_multiqc = harmonise_orfs.out.caller_count_multiqc

    emit:
    // Harmonoised ORF table
    harmonised_orf_table = harmonise_orfs.out.harmonised_orf_table
    // ORFs removed during harmonisation step
    removed_orf_ids = harmonise_orfs.out.removed_orf_ids
    // Annotation statistic files for MultiQC
    annotation_multiqc = orfcaller_multiq.mix(merged_multiqc, caller_count_multiqc).collect()
}
