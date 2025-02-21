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
6. EXPRESSION: Analyze identified ORF expression
7. ANNOTATION: Annotate identified ORFs
*/

//include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SELECTION } from '../subworkflows/SELECTION.nf'
include { ALIGNMENT } from '../subworkflows/ALIGNMENT.nf'
include { RIBOQC } from '../subworkflows/RIBOQC.nf'
include { ORFQUANT } from '../subworkflows/ORFQUANT.nf'
include { PRICE } from '../subworkflows/PRICE.nf'
// Expression and annotation still untested
include { EXPRESSION; EXPRESSION as EXPRESSION_PRICE } from '../subworkflows/EXPRESSION.nf'
//include { ANNOTATION } from '../subworkflows/ANNOTATION.nf'

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

    SELECTION(ch_reads,
              params.bowtie2_index_path,
              params.contaminants_fasta,
              params.keep_sam,
              params.outdir
              )

    rpf_reads = SELECTION.out.rpf_reads

    ALIGNMENT(rpf_reads,
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

    def orfquant_annotation_exists = file(params.orfquant_annotation).isFile()

    RIBOQC(params.orfquant_annotation,
           params.orfquant_annot_package,
           params.package_install_loc,
           params.pandoc_dir,
           orfquant_bams,
           // This was for testing RiboseQC / ORFquant annotation, does not seem to work
           orfquant_annotation_exists,
           params.reference_gtf,
           params.reference_twobit,
           params.outdir,
           params.orfquant_prefix
           // For multiQC and figures
           // contaminants,
           // star_output,
           // samtools_output
           )
    // /*
    if (params.run_orfquant) {
        ORFQUANT(RIBOQC.out.for_orfquant_files,
                 RIBOQC.out.rannot_ch,
                 RIBOQC.out.package_ch,
                 params.package_install_loc,
                 params.pandoc_dir,
                 params.orfquant_prefix,
                 params.outdir
                 )
        
        // This is untested
        // TODO: merge the expression workflows
        ORFQUANT_EXPRESSION = EXPRESSION(ORFQUANT.out.orfquant_orfs,
                                         RIBOQC.out.for_orfquant_files,
                                         params.package_install_loc,
                                         params.outdir
                                         ) 
    }//*/

    if (params.run_price) {
        PRICE(bamlist,
              params.price_index_path,
              params.reference_fasta,
              params.reference_gtf,
              params.outdir,
              params.price_prefix,
              params.gedi_exec_loc
              )
        
         // This is untested
        PRICE_EXPRESSION = EXPRESSION_PRICE(PRICE.out.price_orfs,
                                      RIBOQC.out.for_orfquant_files,
                                      params.package_install_loc,
                                      params.outdir
                                      ) 
    }

/*
    ANNOTATION(
               ORFQUANT.out.orfquant_results,
               PRICE.out.price_results,
               PRICE_EXPRESSION.out.
               ORFQUANT_EXPRESSION.out.
               params.run_orfquant,
               params.run_price
               )
*/

// Run plotting scripts here or somewhere else? Might be better to do it in the seperate processes
// Create multiqc report and add the created plots to it

}
