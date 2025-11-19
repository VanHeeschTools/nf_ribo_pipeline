include { validate_bowtie2_index }                     from "../modules/helperFunctions.nf"
include { trimgalore }                                 from "../modules/trimgalore.nf"
include { fastqc }                                     from "../modules/fastqc.nf"
include { bowtie2; bowtie2_index; contaminants_check } from '../modules/bowtie.nf'

workflow SELECTION {

    // Initial selection of RPF reads for mapping
    // Also outputs QC based on fastq information

    take:
    reads               // Tuple, path to ribo-seq reads with associated sample ID
    bowtie2_index       // Path, precomputed contaminants index for bowtie2
    contaminants_fasta  // Path, fasta file with rRNA, tRNA, and other contaminants
    keep_bam            // Boolean, keep big SAM file for debugging
    outdir              // Path, output directory

    main:
    // Run Trimgalore
    trimgalore(reads, outdir)
    trimmed_reads = trimgalore.out.reads

    // Files for the MultiQC report
    trimgalore_report = trimgalore.out.trimgalore_trimming_report.collect()
    total_reads = trimgalore.out.total_reads.collect()
    removed_reads = trimgalore.out.removed_reads.collect()

    // Validate all bowtie2 index files
    bowtie2_index_check = validate_bowtie2_index(bowtie2_index)

    if (bowtie2_index_check) {
        // If index files already exist, use the provided index prefix
        bowtie2_index_ch = bowtie2_index
        log.info "Using existing Bowtie2 index: ${bowtie2_index_ch}"
    } else {
        log.warn "Some bowtie2 index files are missing. Running bowtie2 indexing."
        // Run bowtie2 - index if none of the index files exist
        bowtie2_index(contaminants_fasta, outdir)
        bowtie2_index_ch = bowtie2_index.out.bowtie2_index_prefix
        log.info "Using created Bowtie2 index: ${bowtie2_index_ch}"
    }

    // Run bowtie2 to filter out contaminants
    bowtie2(bowtie2_index_ch,trimmed_reads, outdir)
    bowtie2_contaminants = bowtie2.out.bowtie_output_files
    rpf_reads = bowtie2.out.filtered_reads

    // Run FASTQC
    fastqc(rpf_reads, outdir)
    fastqc_zip = fastqc.out.fastqc_zip

    // Create QC stats
    contaminants_check(bowtie2_contaminants,
                    keep_bam,
                    outdir)

    // Combine all MultiQC files into one channel
    contaminant_samples = contaminants_check.out.contaminant_samples.collect()
    contaminant_samples_passed = contaminants_check.out.contaminant_samples_passed.collect()
    multiqc_read_samples = trimgalore_report.mix(total_reads,
                                                removed_reads,
                                                fastqc_zip,
                                                contaminant_samples,
                                                contaminant_samples_passed)
                                                .collect()

    emit:
    rpf_reads            // Selected riboseq reads
    multiqc_read_samples // Multiqc input files
}
