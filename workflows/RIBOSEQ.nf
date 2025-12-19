/*
Riboseq Nextflow Pipeline

This pipeline processes ribosome profiling (Riboseq) data, including quality control,
alignment, ORF quantification, and expression analysis.

Main steps:
0.  Validate existence of input samplesheet and input files
1.  SELECTION: Initial read processing and contaminant removal
2.  ALIGNMENT: Align reads to the reference genome
3.  RIBOQC: Quality control of Riboseq data
4.  ORFQUANT: First ORF calling method using ORFquant
5.  PRICE: Second ORF calling method using PRICE
6.  RIBOTIE: Third ORF calling method using RiboTIE
7.  PSITE: Obtain frame of transcripts/ORFs
8.  ANNOTATION: Annotate identified ORFs
9.  EXPRESSION: Analyze identified ORF expression
10. MULTIQC: Create MultiQC report of outputs
*/

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { validateGTF; validate_bowtie2_index; validate_star_index; validate_price_index; copy_samplesheet; collect_output_previous_run } from '../modules/helperFunctions.nf'
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
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files
    // Including if GTF has all required attributes
    if (params.validate_gtf_attributes) {
        validateGTF(params.reference_gtf)
    }

    // Parse samplesheet into channel
    ch_input = channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))
    
    // Filter input samplesheet to keep the fastq files in the rows where the sequence is ribo
    ch_reads = ch_input.filter { meta, _filename_1->
        meta.filetype == "fastq" && meta.sequence == "ribo"
    }

    // Copy input samplesheet to outdir
    copy_samplesheet(params.input, params.outdir)

    // Declare empty multiqc file channel
    multiqc_files = Channel.empty()

    // Run subworkflows

    if (params.run_qc){
        // Quality filtering and contamination removal
        if (params.run_selection){
            SELECTION(
                ch_reads,
                params.bowtie2_index_path,
                params.contaminants_fasta,
                params.keep_bam,
                params.outdir
            )
            rpf_reads = SELECTION.out.rpf_reads
            multiqc_files = multiqc_files.mix(SELECTION.out.multiqc_read_samples)
        } else {
            if (params.run_alignment){
                bowtie_output_files = "${params.outdir}/bowtie2/*/*_filtered.{fastq.gz,fq.gz}"
                rpf_reads = collect_output_previous_run(bowtie_output_files, "sample_id", false, "Bowtie2")
            }
        }

        // Star alignment to reference genome
        if (params.run_alignment){
            ALIGNMENT(
                rpf_reads,
                params.reference_fasta,
                params.star_index_path,
                params.reference_gtf,
                params.outdir
            )
            multiqc_files = multiqc_files.mix(ALIGNMENT.out.star_log_local)

            // Star output channels
            orfquant_bams = ALIGNMENT.out.bam_list
            price_bams = ALIGNMENT.out.bam_list_end2_end.collect()
            ribotie_bams = ALIGNMENT.out.bam_list_end2end_transcriptome
        } else {
            if (params.run_riboseqc){
                // Check if BAM files required for RiboseQC can be found otherwise set to null
                orfquant_bam_files = "${params.outdir}/star/*/*.local.Aligned.sortedByCoord.out.bam"
                orfquant_bams = collect_output_previous_run(orfquant_bam_files, "sample_id", false, "STAR local alignment")
            }
            if (params.run_price){
                // Check if BAM files required for PRICE can be found otherwise set to null
                price_bam_files = "${params.outdir}/star/*/*.end2end.Aligned.sortedByCoord.out.bam"
                price_bams = collect_output_previous_run(price_bam_files, "path", false, "STAR end2end alignment").collect()
            }
            if (params.run_ribotie){
                // Check if BAM files required for RiboTIE can be found otherwise set to null
                ribotie_bam_files ="${params.outdir}/star/*/*.end2end.Aligned.toTranscriptome.out.bam"
                ribotie_bams_search = collect_output_previous_run(ribotie_bam_files, "sample_id", false, "STAR transcriptome end2end alignment")
                ribotie_bams = ribotie_bams_search.map { sample_id, file_list ->
                    [ sample_id, file_list[0] ]
                }
            }
        }

        // RiboseQC run on all local aligned BAM files
        if (params.run_riboseqc){
            html_template = file("${projectDir}/${params.html_template}")
            RIBOQC(
                params.orfquant_annotation,
                params.package_install_loc,
                params.reference_fasta_fai,
                orfquant_bams,
                html_template,
                params.outdir
            )
            for_orfquant_files = RIBOQC.out.for_orfquant_files
            multiqc_files = multiqc_files.mix(RIBOQC.out.multiqc_riboseq)
        } else {
            if (params.run_orfquant || params.run_expression){
                // Check if BAM files required for RiboseQC can be found otherwise set to null
                riboseqc_output_files = "${params.outdir}/riboseqc/*/*_for_ORFquant"
                for_orfquant_files = collect_output_previous_run(riboseqc_output_files, "sample_id", false, "RiboseQC")
            }
        }
    } else{
        if (params.run_orfquant || params.run_expression){
            // Check if BAM files required for RiboseQC can be found otherwise set to null
            riboseqc_output_files = "${params.outdir}/riboseqc/*/*_for_ORFquant"
            for_orfquant_files = collect_output_previous_run(riboseqc_output_files, "sample_id", false, "RiboseQC")
        }
        if (params.run_price){
            // Check if BAM files required for PRICE can be found otherwise set to null
            price_bam_files = "${params.outdir}/star/*/*.end2end.Aligned.sortedByCoord.out.bam"
            price_bams = collect_output_previous_run(price_bam_files, "path", false, "STAR end2end alignment").collect()
        }
        if (params.run_ribotie){
            // Check if BAM files required for RiboTIE can be found otherwise set to null
            ribotie_bam_files ="${params.outdir}/star/*/*.end2end.Aligned.toTranscriptome.out.bam"
            ribotie_bams_search = collect_output_previous_run(ribotie_bam_files, "sample_id", false, "STAR transcriptome end2end alignment")
            ribotie_bams = ribotie_bams_search.map { sample_id, file_list ->
                [ sample_id, file_list[0] ]
            }
        }
    }

    // ORF prediction steps
    if (params.run_orf_prediction) {

        // ORFquant run on merged RiboseQC output
        if (params.run_orfquant){
            ORFQUANT(
                for_orfquant_files,
                params.orfquant_annotation,
                params.reference_gtf,
                params.package_install_loc,
                params.outdir
            )
            orfquant_gtf = ORFQUANT.out.orfquant_orf_gtf
        } else{
            orfquant_output_gtf = "${params.outdir}/orfquant/ORFquant.gtf"
            orfquant_gtf = collect_output_previous_run(orfquant_output_gtf, "path", true, "ORFquant")
        }

        // PRICE run on merged end2end bam files
        if (params.run_price){
            PRICE(
                price_bams,
                params.price_index_path,
                params.reference_fasta,
                params.reference_gtf,
                params.gedi_exec_loc,
                params.outdir
            )          
            price_gtf = PRICE.out.price_orf_gtf
        } else{
            price_output_gtf = "${params.outdir}/price/PRICE.gtf"
            price_gtf = collect_output_previous_run(price_output_gtf, "path", true, "PRICE")
        }

        // Run RiboTIE on transciptome end2end bam files
        if (params.run_ribotie) {
            RIBOTIE(
                ribotie_bams,
                params.reference_fasta,
                params.reference_gtf,
                params.ribotie_min_samples,
                params.outdir
            )
            ribotie_gtf = RIBOTIE.out.ribotie_orf_gtf
        } else {
            ribotie_output_gtf = "${params.outdir}/merged_ribotie/RiboTIE.gtf"
            ribotie_gtf = collect_output_previous_run(ribotie_output_gtf, "path", true, "RIBOTIE.gtf")
        }

        // Combine outputs of ORFcallers into one channel including RiboTIE output
        orfcaller_gtf = price_gtf.mix(orfquant_gtf, ribotie_gtf)

        // Calculate p0 sites in reference CDS and in ORFs
        if (params.run_psite){
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
        } else {
            if (params.run_annotation){
                // Obtain tuple of ORFcaller gtf and ORF - REF psite intersect bed
                search_orfcaller_gtf = "${params.outdir}/annotation/*.gtf"
                if( !file(search_orfcaller_gtf).isEmpty() ) {
                    // Create tuple channel of (gtf, gtf_ref_intersect.bed)
                    orfcaller_gtf_bed = Channel
                        .fromPath(search_orfcaller_gtf)
                        .map { gtf ->
                            def base = gtf.simpleName.replace('.gtf', '')
                            def bed  = file("${gtf.parent}/${base}_ref_intersect.bed")
                            if( !bed.exists() )
                                throw new IllegalStateException("Missing ORF - Ref intersect BED file for ${gtf}")
                            tuple(gtf, bed)
                        }
                        .collect()
                } else {
                    log.info "No ORFcaller GTF files found."
                    orfcaller_gtf_bed = null
                }
                // Search for cds rds file from previous run
                search_ref_cds_rds = "${params.outdir}/annotation/*correct_cds.rds"
                ref_cds_rds = collect_output_previous_run(search_ref_cds_rds, "path", true, "PSITE ref cds rds file")
            }
            // Search for combined ORFcaller psite bed file from previous run
            if (params.run_expression){
                search_orfcaller_psites = "${params.outdir}/annotation/combined_psites.bed"
                orf_gtf_bed = collect_output_previous_run(search_orfcaller_psites, "path", true, "PSITE: combined orfcaller psites bed")
            }
        }

        // Annotate the ORFcaller output gtf files and harmonises them into a single table
        if (params.run_annotation){
            ANNOTATION(
                params.reference_gtf,
                params.package_install_loc,
                orf_gtf_bed,
                ref_cds_rds,
                params.outdir
            )
            harmonised_table = ANNOTATION.out.harmonised_orf_table
            removed_orf_ids = ANNOTATION.out.removed_orf_ids
            multiqc_files = multiqc_files.mix(ANNOTATION.out.annotation_multiqc)
        } else {
            if (params.run_expression){
                harmonised_table_csv = "${params.outdir}/harmonise_orfs/harmonised_orf_table.csv"
                harmonised_table = collect_output_previous_run(harmonised_table_csv, "path", true, "ANNOTATION: harmonised ORF table csv")

                search_removed_orf_ids = "${params.outdir}/harmonise_orfs/removed_orf_ids.txt"
                removed_orf_ids = collect_output_previous_run(search_removed_orf_ids, "path", true, "ANNOTATION: removed ORF ids txt")
            }
        }

        // Calculate expression values for each ORF and add it to the harmonised ORF table
        if (params.run_expression){
            EXPRESSION(
                for_orfquant_files,
                harmonised_table,
                removed_orf_ids,
                orfcaller_psites,
                params.outdir
            )
            multiqc_files = multiqc_files.mix(EXPRESSION.out.multiqc_expression_plot_txt)
        }
    }

    // Creation of MULTIQC report
    if (params.run_multiqc){
        multiqc_config = file("${projectDir}/${params.multiqc_config}")
        MULTIQC(
            multiqc_files.collect(),
            multiqc_config,
            params.outdir
        )
    }
}
