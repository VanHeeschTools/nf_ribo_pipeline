include { riboseqc; sort_bedgraphs; merge_bedgraphs; convert_to_bigwig; create_riboseqc_report } from '../modules/riboseqc.nf'
include { riboseqc_tables; riboseqc_plots } from '../modules/qcplots.nf'

workflow RIBOQC {

    take:
    orfquant_annotation   // Path, ORFquant annotation file
    package_install_loc   // Path, location where BSgenome package is installed
    reference_fasta_fai   // Path, location of reference genome fasta fai file
    orfquant_bams         // List, output from ALIGNMENT subworkflow
    html_template
    outdir                // Path, output directory

    main:
    // 01 - Create riboseqc files
    riboseqc(orfquant_bams,
            outdir,
            orfquant_annotation,
            package_install_loc)

    // 02 - Create p-site tracks
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
                    reference_fasta_fai,
                    outdir)

    // 03 - Create riboseqc tables for MultiQC
    riboseqc_tables(riboseqc.out.riboseqc_all.collect(), outdir)
    riboseqc_inframe_percentages = riboseqc_tables.out.riboseqc_inframe_percentages
    riboseqc_category_counts = riboseqc_tables.out.riboseqc_category_counts

    // 04 - Create periodicity plots
    //riboseqc_plots(riboseqc.out.riboseqc_all.collect(), outdir)


    create_riboseqc_report(
        riboseqc.out.riboseqc_all.collect(),
        html_template,
        outdir
    )

    // Combine into one channel for MultiQC
    multiqc_riboseq = riboseqc_inframe_percentages.mix(riboseqc_category_counts).collect()
                    //riboseqc_plots.out.metagene_plot,
                    //riboseqc_plots.out.periodicity_plot)

    // Obtain ORFquant input files
    // Collect is done in the RIBOSEQ workflow
    for_orfquant_files = riboseqc.out.orfquant_psites

    emit:
    for_orfquant_files     // Files for ORFquant
    multiqc_riboseq        // Files for MultiQC

}
