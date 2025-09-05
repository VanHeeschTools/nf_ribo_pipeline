include { create_template ; parse_genomic_features ; parse_samples ; ribotie_predict_samples ; merge_ribotie_output } from '../modules/ribotie.nf'

workflow RIBOTIE {
    take:
    ribotie_bams        // List, RiboTIE input BAM files
    fasta               // Path, reference FASTA file
    gtf                 // Path, input GTF file
    ribotie_min_samples // Val, min amount of samples the ORF should be found in 
    outdir              // Path, output directory

    main:
    //Create sample template for use in RiboTIE
    create_template(
        ribotie_bams.collect(flat: false),
        fasta,
        gtf,
        outdir
    )
    ribotie_template = create_template.out.template

    // Create h5 database based on genomic features from the gtf 
    parse_genomic_features(
        ribotie_template,
        gtf,
        fasta,
        outdir
    )
    // Define genomic h5 database
    genomic_h5_db = parse_genomic_features.out.h5_path

    // Create h5 database for every sample
    parse_samples(
        genomic_h5_db,
        ribotie_template,
        ribotie_bams,
        gtf,
        fasta,
        outdir
    )
    // Define sample h5 database
    sample_h5 = parse_samples.out.h5_path

    //Run RiboTIE for all samples individually
    ribotie_predict_samples(
        sample_h5,
        genomic_h5_db,
        ribotie_template,
        gtf,
        fasta,
        outdir
    )
    // Define RiboTIE output
    ribotie_orf_csv = ribotie_predict_samples.out.ribotie_orf_csv
    
    // Define RiboTIE generated MultiQC output FOR TESTING
    ribotie_multiqc = ribotie_predict_samples.out.ribotie_multiqc

    // Merge all RiboTIE output files
    merge_ribotie_output(
        ribotie_orf_csv.collect(),
        genomic_h5_db,
        ribotie_min_samples,
        outdir
    )

    // Define RIBOTIE subworkflow output 
    ribotie_orf_gtf = merge_ribotie_output.out.ribotie_merge_gtf
    ribotie_merged = merge_ribotie_output.out.ribotie_merged_csv

    emit:
    ribotie_merged
    ribotie_orf_gtf
    ribotie_multiqc
}
