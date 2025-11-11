include { get_orf_category; harmonise_orfs } from "../modules/annotation.nf"

workflow ANNOTATION {
    take:
    reference_gtf           // Path, input gtf file
    package_install_loc     // Path, Location where BSgenome R package is installed
    orf_gtf_bed             // Path, bed-like file of ORF and ref CDS overlap
    ref_cds_rds             // Path, RDS file of altered reference CDS
    outdir                  // Path, output directory

    main:
    // Load ORFcaller gtf and annotates the ORFs
    get_orf_category(
        orf_gtf_bed,
        reference_gtf,
        ref_cds_rds,
        package_install_loc,
        outdir
    )

    // Collect the annotated ORF tables of all ORFcallers for harmonisation
    annotated_orf_tables = get_orf_category.out.basic_orf_table.collect()

    // Combines the ORFcaller annotated csv files into one harmonised ORF table
    harmonise_orfs(
        annotated_orf_tables,
        outdir
    )

    // Define subworkflow output
    harmonised_orf_table = harmonise_orfs.out.harmonised_orf_table
    removed_orf_ids = harmonise_orfs.out.removed_orf_ids

    // Annotation multiqc output files
    orfcaller_multiq = harmonise_orfs.out.orfcaller_multiq
    merged_multiqc = harmonise_orfs.out.merged_multiqc
    caller_count_multiqc = harmonise_orfs.out.caller_count_multiqc
    
    annotation_multiqc = orfcaller_multiq.mix(merged_multiqc, caller_count_multiqc).collect()

    emit:
    harmonised_orf_table
    removed_orf_ids
    annotation_multiqc
}
