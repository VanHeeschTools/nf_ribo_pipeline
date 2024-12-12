include { trimming } from "../modules/fastp.nf"
include { fastqc } from "../modules/fastqc.nf"
include { bowtie2; bowtie2_index } from '../modules/bowtie.nf'
include { contaminants_check } from '../modules/samtools.nf'

workflow SELECTION {

    // Initial selection of RPF reads for mapping
    // Also outputs QC based on fastq information

    take:
    reads               // Riboseq reads with associated sample ID in tuple format
    bowtie2_index       // Precomputed contaminants index for bowtie2
    contaminants_fasta  // Fasta file with rRNA, tRNA, and other contaminants
    keep_sam            // Boolean, keep big SAM file for debugging
    outdir              // Output directory

    main:
    // Run FASTP
    trimming(reads, outdir)

    // Run FASTQC
    fastqc(trimming.out.reads, outdir)

    // Define the Bowtie2 index file extensions
    def bowtie2_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

    // Check if any of the Bowtie2 index files exist
    def index_files_exist = bowtie2_extensions.every { ext -> file(bowtie2_index + ext).exists() }
    log.info "All bowtie2 index files exist: ${index_files_exist}"

    if (!index_files_exist) {
        log.warn "Some bowtie2 index files are missing. Running bowtie2 indexing."
        // Run bowtie2 - index if none of the index files exist
        bowtie2_index(contaminants_fasta, outdir)
        bowtie2_index_ch = bowtie2_index.out.bowtie2_index
    } else {
        // If index files already exist, use the provided index prefix
        bowtie2_index_ch = bowtie2_index
        log.info "Using existing Bowtie index: ${bowtie2_index_ch}"
    }

    // Run bowtie2
    bowtie2(bowtie2_index_ch,trimming.out.reads, outdir)

    // Create QC stats
    contaminants_check(trimming.out.reads,
                       bowtie2.out.filtered_reads,
                       bowtie2.out.sam_file,
                       bowtie2_index,
                       outdir,
                       keep_sam)

    // Gather output tuples into single list
    rpf_reads = bowtie2.out.filtered_reads
    contaminants_bowtie2 = contaminants_check.out.contaminants_bowtie2

    emit:
    rpf_reads            // Selected riboseq reads
    contaminants_bowtie2 // Contaminant stats per sample
}
