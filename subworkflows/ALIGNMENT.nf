include { star_index ; star_local ; star_end_to_end } from '../modules/star.nf'
include { samtools ; samtools as samtools_end2end } from '../modules/samtools.nf'

workflow ALIGNMENT {
    take:
    rpf_reads          // Path, output from SELECTION subworkflow
    genome             // Path, reference genome used for STAR index
    star_index_path    // Path, location of precomputed STAR index
    gtf                // Path, reference gtf file
    run_orf_prediction // Bool, True if orf callers should be run
    outdir             // Path, output directory

    main:
    // List of files for which to check existence
    def star_index_files = [
        'chrLength.txt',
        'chrStart.txt',
        'geneInfo.tab',
        'SA',
        'sjdbList.fromGTF.out.tab',
        'chrNameLength.txt',
        'exonGeTrInfo.tab',
        'Genome',
        'SAindex',
        'sjdbList.out.tab',
        'chrName.txt',
        'exonInfo.tab',
        'genomeParameters.txt',
        'sjdbInfo.txt',
        'transcriptInfo.tab',
    ]

    // Check if the STAR index files exist
    def star_index_files_exist = star_index_files.every { filename ->
        file("${star_index_path}/${filename}").exists()
    }
    log.info("All STAR index files exist: ${star_index_files_exist}")

    if (star_index_files_exist == false) {
        log.warn("Some STAR index files are missing. Running STAR indexing.")
        // Run STAR - index if any of the index files are missing
        star_index(genome, gtf, outdir)
        star_index_ch = star_index.out.star_index_path
    } else {
        star_index_ch = "${star_index_path}"
        log.info("Using existing STAR index: ${star_index_path}")
    }

    // Run STAR local mode
    star_local(
        rpf_reads,
        outdir,
        gtf,
        star_index_ch,
    )

    // Sort output BAM file 
    samtools(star_local.out.bams, outdir)
    bam_list = samtools.out.sorted_bam
    star_log_local = star_local.out.star_log_local

    // Only run the end2end STAR run if the ORF calling will be run
    // This is because at the moment only the ORFquant BAM file is used for MultiQC output
    if (run_orf_prediction) {
        // Run STAR end2end mode
        star_end_to_end(
            rpf_reads,
            outdir,
            gtf,
            star_index_ch,
        )

        bam_list_end2end_transcriptome = star_end_to_end.out.bams_end2end_transcriptome

        // Sort output BAM file 
        samtools_end2end(star_end_to_end.out.bams_end2end, outdir)
        star_log_end_to_end = star_end_to_end.out.star_log_end_to_end

        // Obtain all STAR end2end sorted BAM file paths and change path to string
        price_paths = samtools_end2end.out.bam_files
            .collect()
            .flatten()
            .map { it -> it.toString() }

        //TODO: Remove this part and use price_paths as input for the merge_price_bams step
        // Store BAM file paths to workdir
        price_filelist = price_paths.collectFile(
            name: 'bams.bamlist',
            newLine: true,
            sort: true,
        )
    } else {
        star_log_end_to_end = null
        bam_list_end2end_transcriptome = null
        price_filelist = null
    }

    emit:
    star_log_local                 // star output log file for local run
    star_log_end_to_end            // star output log file for end2end run
    bam_list                       // bam files for ORFquant
    bam_list_end2end_transcriptome // bam list for RiboTIE
    price_filelist                 // list for PRICE
}
