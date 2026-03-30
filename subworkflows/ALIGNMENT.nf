include { validate_star_index }                       from "../modules/helperFunctions.nf"
include { star_index ; star_local ; star_end_to_end } from '../modules/star.nf'
include { samtools ; samtools as samtools_end2end; samtools as samtools_transcriptome  }   from '../modules/samtools.nf'

workflow ALIGNMENT {
    take:
    rpf_reads          // Path, output from SELECTION subworkflow
    genome             // Path, reference genome used for STAR index
    star_index_path    // Path, location of precomputed STAR index
    gtf                // Path, reference gtf file

    main:

    // Validate all STAR index files
    star_index_check =  validate_star_index(params.star_index_path)

    // Create STAR index if any of the index files is missing
    if (star_index_check) {
        star_index_ch = "${star_index_path}"
        log.info("Using existing STAR index: ${star_index_ch}")
    } else {
        log.warn("Some STAR index files are missing. Running STAR indexing.")
        // Run STAR - index if any of the index files are missing
        star_index(genome, gtf)
        star_index_ch = star_index.out.star_index_path
    }

    // Run STAR local mode
    star_local(
        rpf_reads,
        gtf,
        star_index_ch,
    )

    // Sort local BAM file 
    samtools(star_local.out.bams)
    bam_list = samtools.out.sorted_bam
    star_log_local = star_local.out.star_log_local

    // Run STAR end2end mode
    star_end_to_end(
        rpf_reads,
        gtf,
        star_index_ch,
    )

    bam_list_end2end_transcriptome = star_end_to_end.out.bams_end2end_transcriptome

    // Sort end2end BAM file 
    samtools_end2end(star_end_to_end.out.bams_end2end)
    star_log_end_to_end = star_end_to_end.out.star_log_end_to_end

    // Sort end2end transcriptome BAM file 
    samtools_transcriptome(bam_list_end2end_transcriptome)
    bam_list_end2end_transcriptome_sorted = samtools_transcriptome.out.sorted_bam

    // Obtain all STAR end2end sorted BAM file paths and change path to string
    bam_list_end2_end = samtools_end2end.out.bam_files
        .collect()
        .flatten()
        .map { it -> it.toString() }

    emit:
    star_log_local                 // star output log file for local run
    star_log_end_to_end            // star output log file for end2end run
    bam_list                       // bam files for ORFquant
    bam_list_end2_end              // bam files for PRICE
    bam_list_end2end_transcriptome_sorted // bam list for RiboTIE
}
