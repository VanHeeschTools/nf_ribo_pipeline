/*
Riboseq Nextflow Pipeline

This pipeline processes ribosome profiling (Riboseq) data, including quality control,
alignment, ORF quantification, and expression analysis.

Main steps:
1. SELECTION: Initial read processing and contaminant removal
2. ALIGNMENT: Align reads to the reference genome
3. RIBOQC: Quality control of Riboseq data
4. ORFQUANT: ORF quantification
5. PRICE: Alternative ORF calling method
6. PSITE: Obtain frame of transcripts/ORFs
7. ANNOTATION: Annotate identified ORFs
8. EXPRESSION: Analyze identified ORF expression

*/

//include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { SELECTION }  from '../subworkflows/SELECTION.nf'
include { ALIGNMENT }  from '../subworkflows/ALIGNMENT.nf'
include { RIBOQC }     from '../subworkflows/RIBOQC.nf'
include { ORFQUANT }   from '../subworkflows/ORFQUANT.nf'
include { PRICE }      from '../subworkflows/PRICE.nf'
include { PSITE }      from '../subworkflows/PSITE.nf'
include { ANNOTATION } from '../subworkflows/ANNOTATION.nf'
include { EXPRESSION } from '../subworkflows/EXPRESSION.nf'
include { MULTIQC }    from '../modules/multiqc.nf'

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
    // Including if GTF has all required attributes

    // TODO: handle annotations from samplesheet

    // Set multiqc file channel
    multiqc_files = Channel.empty()

    SELECTION(ch_reads,
              params.bowtie2_index_path,
              params.contaminants_fasta,
              params.keep_sam,
              params.outdir
              )
    rpf_reads = SELECTION.out.rpf_reads
    multiqc_files = multiqc_files.mix(SELECTION.out.multiqc_read_samples)

    ALIGNMENT(rpf_reads,
              params.reference_fasta,
              params.star_index_path,
              params.reference_gtf,
              params.outdir
              )
    multiqc_files = multiqc_files.mix(ALIGNMENT.out.star_log_local)

    orfquant_bams = ALIGNMENT.out.bam_list
    //price_bams = ALIGNMENT.out.bam_list_end2end // Currently unused
    bamlist = ALIGNMENT.out.price_filelist

    def orfquant_annotation_exists = file(params.orfquant_annotation).isFile()
    //def orfquant_annotation_exists = false

    RIBOQC(params.orfquant_annotation,
           params.orfquant_annot_package,
           params.package_install_loc,
           params.pandoc_dir,
           orfquant_bams,

           // This was for testing RiboseQC / ORFquant annotation, does not seem to work
           orfquant_annotation_exists,
           params.reference_gtf,
           params.reference_twobit,
           params.reference_fasta,
           params.outdir
           )
    multiqc_files = multiqc_files.mix(RIBOQC.out.multiqc_riboseq)

    // ORF prediction steps
    if (params.run_orf_prediction ) {
        ORFQUANT(RIBOQC.out.for_orfquant_files,
                RIBOQC.out.rannot_ch,
                RIBOQC.out.package_ch,
                params.package_install_loc,
                params.pandoc_dir,
                params.outdir
                )

        PRICE(bamlist,
            params.price_index_path,
            params.reference_fasta,
            params.reference_gtf,
            params.gedi_exec_loc,
            params.outdir
            )

        // Combine outputs of ORFcallers into one channel
        orfcaller_gtf = PRICE.out.price_orf_gtf.mix(ORFQUANT.out.orfquant_orf_gtf)
        orfcaller_output = PRICE.out.price_orfs.mix(ORFQUANT.out.orfquant_orfs)

        PSITE(orfcaller_gtf,
            params.reference_gtf,
            RIBOQC.out.for_orfquant_files,
            params.outdir
            ) 

        orfcaller_psites = PSITE.out.orfcaller_psites
        ref_psites = PSITE.out.ref_psites


        ANNOTATION(orfcaller_output,
                orfcaller_psites,
                ref_psites,
                params.reference_gtf,
                params.package_install_loc,
                params.orfquant_annot_package,
                params.outdir
                )

        EXPRESSION(RIBOQC.out.for_orfquant_files,
                ANNOTATION.out.harmonised_orf_table,
                ANNOTATION.out.removed_orf_ids,
                orfcaller_psites,
                params.outdir)

    }
    multiqc_config = file("${projectDir}/${params.multiqc_config}")
    MULTIQC(multiqc_files.collect(),
            multiqc_config,
            params.outdir)

// Run plotting scripts here or somewhere else? Might be better to do it in the seperate processes
// Create multiqc report and add the created plots to it

}
