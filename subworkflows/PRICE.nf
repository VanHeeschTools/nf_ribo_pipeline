include { price; price_index; price_to_gtf; merge_price_bams } from '../modules/price.nf'

workflow PRICE {

    take:
    bamlist           // List of PRICE input BAM files
    price_index_path  // Path to the PRICE index file
    fasta             // Reference FASTA file
    gtf               // Input GTF file
    gedi_exec_loc     // Location of local GEDI installation
    outdir            // Output directory

    main:
    // Check if PRICE index exists
    def price_index_exists = file("${price_index_path}/PRICE_index.oml").exists()
    log.info "PRICE index exists: ${price_index_exists}"

    if (!price_index_exists) {
        // Create PRICE annotation
        log.warn "PRICE index file missing. Running PRICE indexing."
        price_index(fasta,
                    gtf,
                    gedi_exec_loc)
        price_index_ch = price_index.out.price_index

    } else {
        price_index_ch = Channel.value(file("${price_index_path}/PRICE_index.oml"))
        log.info "Using existing PRICE index: ${price_index_ch}"
    }
    
    // Merge bam files before running PRICE
    merge_price_bams(bamlist,
                     outdir)

    // Run PRICE
    price(merge_price_bams.out.merged_end2end_bam,
      price_index_ch,
      gedi_exec_loc,
      outdir)
    
    // Turn PRICE output into gtf format
    price_to_gtf(
      price.out.price_orfs,
      outdir
    )

    // Set output variables
    price_orfs = price.out.price_orfs
    price_orf_gtf = price_to_gtf.out.price_gtf

    emit:
    price_orfs
    price_orf_gtf

}
