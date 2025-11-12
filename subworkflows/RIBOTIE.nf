include { create_template ; parse_genomic_features ; parse_samples ; ribotie_predict_samples ; merge_ribotie_output; ribotie_add_stop } from '../modules/ribotie.nf'

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
        gtf,
        fasta,
        outdir
    )
    ribotie_template = create_template.out.template

    //TODO: see if this step can be run without the template so the process can start earlier
    // Create h5 database based on genomic features from the gtf 
    parse_genomic_features(
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
    
    // Merge all RiboTIE output files
    merge_ribotie_output(
        ribotie_orf_csv.collect(),
        genomic_h5_db,
        ribotie_min_samples,
        outdir
    )

    // Define RIBOTIE subworkflow output 
    ribotie_gtf_no_stop = merge_ribotie_output.out.ribotie_merged_gtf

    // Add stop codon coords to predicted ORFs exon boundary aware
    ribotie_add_stop(
        ribotie_gtf_no_stop,
        gtf,
        outdir
    )

    ribotie_orf_gtf = ribotie_add_stop.out.ribotie_gtf

    emit:
    ribotie_orf_gtf
}
