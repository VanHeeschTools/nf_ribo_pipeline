//include { ref_psites; sample_psites} from "../modules/process_bed.nf"
include { annotate_orfs; annotate_orfs as annotate_price_orfs; harmonise_orfs} from "../modules/annotation.nf"

workflow ANNOTATION {

    take:
    orfquant_results
    orfquant_ppm_matrix
    price_results
    price_ppm_matrix
    reference_gtf
    run_orfquant
    run_price
    annotation_provider
    gencode_uniprot_file
    uniprot_protein_fasta_loc
    package_install_loc
    orfquant_annot_package
   

    main:
    if (run_orfquant){
        orfcaller = "orfquant"
        annotate_orfs(
            orfquant_results,
            reference_gtf,
            orfcaller,
            annotation_provider,
            gencode_uniprot_file,
            uniprot_protein_fasta_loc,
            package_install_loc,
            orfquant_annot_package
        )
    }
    
    if (run_price){
        orfcaller = "price"
        annotate_price_orfs(
            price_results,
            reference_gtf,
            orfcaller,
            annotation_provider,
            gencode_uniprot_file,
            uniprot_protein_fasta_loc,
            package_install_loc,
            orfquant_annot_package
        )
    }

    //if (run_orfquant && run_price){
    //    harmonise_orfs
    //}



    //emit:
    //sqharmonised_orf_table


}
