include { riboseqc; sort_bedgraphs; merge_bedgraphs; convert_to_bigwig; create_riboseqc_report; riboseqc_tables } from '../modules/riboseqc.nf'
workflow RIBOQC {

    take:
    orfquant_annotation      // Path, ORFquant annotation file
    package_install_loc      // Path, location where BSgenome package is installed
    readlength_choice_method // Val, filter choice of RiboseQC should be "max_coverage" or "all"
    reference_fasta_fai      // Path, location of reference genome fasta fai file
    orfquant_bams            // List, output from ALIGNMENT subworkflow
    html_template            // Path, location of RiboseQC html report template

    main:
    // Create RiboseQC files
    riboseqc(orfquant_bams,
            orfquant_annotation,
            package_install_loc,
            readlength_choice_method)

    // Create p-site tracks
    // Sort each bedgraph file
    riboseqc_bedgraphs = riboseqc.out.bedgraphs.flatten()
    sort_bedgraphs(riboseqc_bedgraphs)
    sorted_ch = sort_bedgraphs.out.sorted_bedgraph.collect()

    // Merge bedgraphs into groups
    merge_bedgraphs(sorted_ch)
    merged_bedgraphs = merge_bedgraphs.out.flatten()

    // Mix in unmerged files to create seperate bigwig files
    list_of_bedgraphs = merged_bedgraphs.mix(sorted_ch).flatten()
    // Convert bedgraph to bigwig
    convert_to_bigwig(list_of_bedgraphs,
                    reference_fasta_fai)

    // Create RiboseQC tables for MultiQC
    riboseqc_tables(riboseqc.out.riboseqc_all.collect())
    riboseqc_inframe_percentages = riboseqc_tables.out.riboseqc_inframe_percentages
    riboseqc_category_counts = riboseqc_tables.out.riboseqc_category_counts

    // Creat RiboseQC output html report
    create_riboseqc_report(
        riboseqc.out.riboseqc_all.collect(),
        html_template
    )

    emit:
    // Obtain ORFquant input files, collect is used in the RIBOSEQ workflow
    for_orfquant_files = riboseqc.out.orfquant_psites
    // Combine into one channel for MultiQC
    multiqc_riboseq = riboseqc_inframe_percentages.mix(riboseqc_category_counts).collect()

}
