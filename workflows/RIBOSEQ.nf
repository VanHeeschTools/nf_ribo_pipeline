/*
Riboseq Nextflow Pipeline

This pipeline processes ribosome profiling (Riboseq) data, including quality control,
alignment, ORF quantification, and expression analysis.

Main steps:
1.  SELECTION: Initial read processing and contaminant removal
2.  ALIGNMENT: Align reads to the reference genome
3.  RIBOQC: Quality control of Riboseq data
4.  ORFQUANT: ORF quantification
5.  PRICE: Alternative ORF calling method
6.  RiboTIE: Third ORF calling method
7.  PSITE: Obtain frame of transcripts/ORFs
8.  ANNOTATION: Annotate identified ORFs
9.  EXPRESSION: Analyze identified ORF expression
10. MULTIQC: Create MultiQC report of outputs
*/

//include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { validateGTF } from '../modules/helperFunctions.nf'
include { SELECTION   } from '../subworkflows/SELECTION.nf'
include { ALIGNMENT   } from '../subworkflows/ALIGNMENT.nf'
include { RIBOQC      } from '../subworkflows/RIBOQC.nf'
include { ORFQUANT    } from '../subworkflows/ORFQUANT.nf'
include { PRICE       } from '../subworkflows/PRICE.nf'
include { RIBOTIE     } from '../subworkflows/RIBOTIE.nf'
include { PSITE       } from '../subworkflows/PSITE.nf'
include { ANNOTATION  } from '../subworkflows/ANNOTATION.nf'
include { EXPRESSION  } from '../subworkflows/EXPRESSION.nf'
include { MULTIQC     } from '../modules/multiqc.nf'

workflow RIBOSEQ {

    // Validate input parameters
    // validateParameters()

    // Print summary of supplied parameters
    // log.info paramsSummaryLog(workflow)

    // Validate files
    // Including if GTF has all required attributes
    if (params.validate_gtf_attributes) {
        validateGTF(params.reference_gtf)
    }
    // TODO: handle annotations from samplesheet

    // Define input channel
    ch_input = Channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.sample_id = row.sample_id
            meta.subject_id = row.subject_id
            meta.sample_type = row.sample_type
            meta.sequence_type = row.sequence_type
            meta.file_type = row.file_type
            // Use 'realpath' to resolve potential symlinks
            def resolvedPath = "realpath ${row.filename_1}".execute().text.trim()
            [meta, file(resolvedPath)]
        }

    // Filter for ribo-seq fastq files
    ch_reads = ch_input.filter { meta, _fastq ->
        meta.file_type == "fastq" && meta.sequence_type == "ribo"
    }

    // Declare empty multiqc file channel
    multiqc_files = Channel.empty()

    // Quality filtering and contamination removal
    SELECTION(
        ch_reads,
        params.bowtie2_index_path,
        params.contaminants_fasta,
        params.keep_bam,
        params.outdir
    )
    rpf_reads = SELECTION.out.rpf_reads
    multiqc_files = multiqc_files.mix(SELECTION.out.multiqc_read_samples)

    // Star alignment to reference genome
    ALIGNMENT(
        rpf_reads,
        params.reference_fasta,
        params.star_index_path,
        params.reference_gtf,
        params.run_orf_prediction,
        params.outdir
    )
    multiqc_files = multiqc_files.mix(ALIGNMENT.out.star_log_local)

    // Star output channels
    orfquant_bams = ALIGNMENT.out.bam_list
    bamlist = ALIGNMENT.out.price_filelist
    ribotie_bams = ALIGNMENT.out.bam_list_end2end_transcriptome

    // RiboseQC run on all local aligned BAM files
    RIBOQC(
        params.orfquant_annotation,
        params.package_install_loc,
        params.reference_fasta_fai,
        orfquant_bams,
        params.outdir
    )
    multiqc_files = multiqc_files.mix(RIBOQC.out.multiqc_riboseq)

    // ORF prediction steps
    if (params.run_orf_prediction) {
        // ORFquant run on merged RiboseQC output
        ORFQUANT(
            RIBOQC.out.for_orfquant_files,
            params.orfquant_annotation,
            params.reference_gtf,
            params.package_install_loc,
            params.outdir
        )

        // PRICE run on merged end2end bam files
        PRICE(
            bamlist,
            params.price_index_path,
            params.reference_fasta,
            params.reference_gtf,
            params.gedi_exec_loc,
            params.outdir
        )

        // Run RiboTIE only if run_ribotie parameter is set to TRUE
        if (params.run_ribotie) {
            // Run RiboTIE on transciptome end2end bam files
            RIBOTIE(
                ribotie_bams,
                params.reference_fasta,
                params.reference_gtf,
                params.ribotie_min_samples,
                params.outdir
            )

            // Combine outputs of ORFcallers into one channel including RiboTIE output
            orfcaller_gtf = PRICE.out.price_orf_gtf.mix(ORFQUANT.out.orfquant_orf_gtf, RIBOTIE.out.ribotie_orf_gtf)
        }
        else {
            // Combine outputs of ORFcallers into one channel excluding RiboTIE output
            orfcaller_gtf = PRICE.out.price_orf_gtf.mix(ORFQUANT.out.orfquant_orf_gtf)
        }
        
        // Calculate p0 sites in reference CDS and in ORFs
        PSITE(
            orfcaller_gtf,
            params.reference_gtf,
            params.reference_protein_fa,
            params.package_install_loc,
            params.outdir
        )

        // Merged ORFcallers p0 psites, expression input
        orfcaller_psites = PSITE.out.orfcaller_psites

        // ORFcaller gtf file plus reference p-site overlap bed file, annotation input
        orf_gtf_bed = PSITE.out.orf_gtf_bed
        
        // Altered reference cds rds file, annotation input
        ref_cds_rds = PSITE.out.ref_cds_rds

        // Annotate the ORFcaller output gtf files and harmonises them into a single table
        ANNOTATION(
            params.reference_gtf,
            params.package_install_loc,
            orf_gtf_bed,
            ref_cds_rds,
            params.outdir
        )
        multiqc_files = multiqc_files.mix(ANNOTATION.out.annotation_multiqc)

        // Calculate expression values for each ORF and add it to the harmonised ORF table
        EXPRESSION(
            RIBOQC.out.for_orfquant_files,
            ANNOTATION.out.harmonised_orf_table,
            ANNOTATION.out.removed_orf_ids,
            orfcaller_psites,
            params.outdir
        )
        multiqc_files = multiqc_files.mix(EXPRESSION.out.multiqc_expression_plot_txt)
    }

    // Creation of MULTIQC report
    multiqc_config = file("${projectDir}/${params.multiqc_config}")
    MULTIQC(
        multiqc_files.collect(),
        multiqc_config,
        params.outdir
    )
}
