/*
Riboseq Nextflow Pipeline

This pipeline processes ribosome profiling (Riboseq) data, including quality control,
alignment, ORF quantification, and expression analysis.

Main steps:
1. SELECTION: Initial read processing and contaminant removal
2. ALIGNMENT: Align reads to the reference genome
3. RIBOQC: Quality control of Riboseq data
4. ORFQUANT: ORF quantification (optional)
5. PRICE: Alternative ORF calling method (optional)
6. ANNOTATION: Annotate identified ORFs
7. EXPRESSION: Analyze ORF expression
*/

//include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SELECTION } from '../subworkflows/selection.nf'
include { ALIGNMENT } from '../subworkflows/alignment.nf'
//include { RIBOQC } from '../subworkflows/riboqc.nf'
//include { ORFQUANT } from '../subworkflows/orfquant.nf'
include { PRICE } from '../subworkflows/price.nf'
/* include { EXPRESSION } from '../subworkflows/expression.nf'
include { ANNOTATION } from '../subworkflows/annotation.nf' */

// Define input channel
ch_input = Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row -> 
        def meta = [:]
        meta.sample_id = row.sample_id
        meta.subject_id = row.subject_id
        meta.sample_type = row.sample_type
        meta.sequence_type = row.sequence_type
        meta.file_type = row.file_type
        // Use 'realpath' to resolve the symlink
        def resolvedPath = "realpath ${row.filename_1}".execute().text.trim()
        
        [ meta, file(resolvedPath) ]
    }

// Filter for ribo-seq fastq files
ch_reads = ch_input.filter { meta, fastq ->
        meta.file_type == "fastq" && meta.sequence_type == "ribo"
    }

workflow RIBOSEQ {

    // Validate input parameters
    // validateParameters()

    // Print summary of supplied parameters
    // log.info paramsSummaryLog(workflow)

    // Validate files

    // TODO: handle annotations from samplesheet

    SELECTION(
        ch_reads,
        params.bowtie2_index_path,
        params.contaminants_fasta,
        params.keep_sam,
        params.outdir
        )

    rpf_reads = SELECTION.out.rpf_reads

    ALIGNMENT(
        rpf_reads,
        params.reference_fasta,
        params.star_index_path,
        params.reference_gtf,
        params.run_price,
        params.run_orfquant,
        params.outdir
        )

    orfquant_bams = ALIGNMENT.out.bam_list
    price_bams = ALIGNMENT.out.bam_list_end2end
    bamlist = ALIGNMENT.out.price_filelist

    //def orfquant_annotation_exists = params.orfquant_annotation.exists()
/*
    RIBOQC(
           orfquant_annotation_exists,
           orfquant_bams,
           params.reference_gtf,
           params.reference_twobit,
           params.outdir,
           params.orfquant_prefix,
           contaminants,
           star_output,
           samtools_output
           )


    if (params.run_orfquant) {
        ORFQUANT(RIBOQC.out.for_orfquant_files)
    }
*/
    if (params.run_price) {
        PRICE(
            bamlist,
            params.price_index_path,
            params.reference_fasta,
            params.reference_gtf,
            params.outdir,
            params.price_prefix,
            params.gedi_exec_loc
            )
    }
/*
    ANNOTATION(
        ORFQUANT.out.orfquant_results,
        PRICE.out.price_results,
        gtf_ch,
        params.run_orfquant,
        params.run_price
    )
    
    EXPRESSION(
        ANNOTATION.out.annotated_orfs,
        ALIGNMENT.out.bam_files,
        params.run_orfquant,
        params.run_price
    ) */

}
