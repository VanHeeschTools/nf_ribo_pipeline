include { price; price_index; price_to_gtf } from '../modules/price.nf'

workflow PRICE {

    take:
    bamlist
    price_index_path
    fasta
    gtf               // Transc
    outdir            // Output directory
    price_prefix      // PRICE name
    gedi_exec_loc     // Location of local GEDI installation

    main:
    // Check if PRICE index exists
    def price_index_exists = file("${price_index_path}/${price_prefix}.oml").exists()
    log.info "PRICE index exists: ${price_index_exists}"

    if (!price_index_exists) {
        // Create PRICE annotation
        log.warn "PRICE index file missing. Running PRICE indexing."
        price_index(fasta,
                    gtf,
                    price_prefix,
                    gedi_exec_loc)
        price_index_ch = price_index.out.price_index

    } else {
        price_index_ch = Channel.value(file("${price_index_path}/${price_prefix}.oml"))
        log.info "Using existing PRICE index: ${price_index_ch}"
    }
    
    price(bamlist,
      price_index_ch,
      price_prefix,
      gedi_exec_loc,
      outdir)
    
    price_to_gtf(
      price.out.price_orfs,
      outdir
    )

    price_orfs = price.out.price_orfs
    price_orf_gtf = price_to_gtf.out.price_gtf

    emit:
    price_orfs
    price_orf_gtf

}
