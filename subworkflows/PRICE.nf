include { validate_price_index }                               from "../modules/helperFunctions.nf"
include { price; price_index; price_to_gtf; merge_price_bams } from '../modules/price.nf'

// Merge end2end BAM files and run PRICE
workflow PRICE {

  take:
  bamlist           // List, PRICE input BAM files
  price_index       // Path, PRICE index file
  fasta             // Path, reference FASTA file
  gtf               // Path, input GTF file
  gedi_exec_loc     // Path, location of local GEDI installation
  outdir            // Path, output directory

  main:

  // Validate price index file
  price_index_check = validate_price_index(price_index)

  // Create PRICE index if it is not found
  if (price_index_check) {
      price_index_ch = Channel.value(file("${price_index}/PRICE_index.oml"))
      log.info "Using existing PRICE index: ${price_index_ch}"
  } else {
      // Create PRICE annotation
      log.warn "PRICE index file missing. Running PRICE indexing."
      price_index(fasta,
                  gtf,
                  gedi_exec_loc)
      price_index_ch = price_index.out.price_index
  }
  
  // Merge bam files before running PRICE
  merge_price_bams(bamlist,
    outdir)

  // Run PRICE
  price(merge_price_bams.out.merged_end2end_bam,
    price_index_ch,
    gedi_exec_loc,
    outdir
  )
  
  // Turn PRICE output into gtf format
  price_to_gtf(
    price.out.price_orfs,
    outdir
  )

  // Define PRICE subworkflow output
  price_orfs = price.out.price_orfs
  price_orf_gtf = price_to_gtf.out.price_gtf

  emit:
  price_orfs
  price_orf_gtf

}
