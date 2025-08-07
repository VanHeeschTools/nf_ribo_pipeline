include { create_template; parse_genomic_features; parse_samples; ribotie_predict_samples; merge_ribotie_output } from '../modules/ribotie.nf'

workflow RIBOTIE{

take:

ribotie_bams        // List of PRICE input BAM files
fasta               // Reference FASTA file
gtf                 // Input GTF file
ribotie_min_samples // Min amount of samples the ORF should be found in 
outdir              // Output directory

main:

//test = ribotie_bams.collect(flat: false)
//test.view()

//Create sample template for use in RiboTIE
create_template(ribotie_bams.collect(flat: false), 
                fasta,
                gtf,
                outdir)
ribotie_template = create_template.out.template

// Create h5 database based on genomic features from the gtf 
parse_genomic_features(ribotie_template,
                       gtf, 
                       fasta,
                       outdir)

genomic_h5_db = parse_genomic_features.out.h5_path

// Create h5 database for every sample
parse_samples(genomic_h5_db,
              ribotie_template,
              ribotie_bams,
              gtf, 
              fasta,
              outdir)
sample_h5 = parse_samples.out.h5_path
    
//Run RiboTIE for all samples individually
ribotie_predict_samples(sample_h5,
                       genomic_h5_db,
                       ribotie_template,
                       gtf, 
                       fasta,
                       outdir
)
//ribotie_orf_gtf = ribotie_predict_samples.out.ribotie_orf_gtf
ribotie_orf_csv = ribotie_predict_samples.out.ribotie_orf_csv

ribotie_multiqc = ribotie_predict_samples.out.ribotie_multiqc


//merge_ribotie_output(ribotie_orf_csv.collect(), 
//                     ribotie_orf_gtf.collect())

merge_ribotie_output(ribotie_orf_csv.collect(),
                    genomic_h5_db,
                    ribotie_min_samples)

ribotie_orf_gtf = merge_ribotie_output.out.ribotie_merge_gtf
ribotie_merged = merge_ribotie_output.out.ribotie_merged_csv


emit:
ribotie_merged
ribotie_orf_gtf
ribotie_multiqc

}