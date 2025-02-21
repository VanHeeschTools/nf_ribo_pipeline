//include { ref_psites; sample_psites} from "../modules/process_bed.nf"
include { annotate_orfs} from "../modules/annotation.nf"

workflow ANNOTATION {

    take:
    orfquant_results
    price_results
    PRICE_EXPRESSION.out
    ORFQUANT_EXPRESSION.out
    run_orfquant
    run_price
   

    main:

    annotate_orfs()


    emit:
    harmonised_orf_table


}
