include { annotate_orfs; harmonise_orfs } from "../modules/annotation.nf"

workflow ANNOTATION {

    take:
    orfcaller_output
    orfcaller_psites
    ref_psites
    reference_gtf
    package_install_loc
    orfquant_annot_package
    run_ribotie
    outdir
   
    main:
    // Parses and annotates the output of the ORFcallers
    annotate_orfs(
        orfcaller_output,
        orfcaller_psites,
        ref_psites,
        reference_gtf,
        package_install_loc,
        orfquant_annot_package,
        outdir
    )

     // Collect ORF tables and sort them alphabetically
    if (run_ribotie){
        harmonise_input = annotate_orfs.out.basic_orf_table
        .collect()
        .map { files ->
            def sorted = files.sort { it.name }
            tuple( sorted[0], sorted[1], sorted[2] )
        }
    } else {
        harmonise_input = annotate_orfs.out.basic_orf_table
        .collect()
        .map { files ->
            def sorted = files.sort { it.name }
            tuple( sorted[0], sorted[1], null )
        }
    }

    // Combines the ORFcaller annotated csv files
    harmonise_orfs(harmonise_input,
                   outdir)

    harmonised_orf_table = harmonise_orfs.out.harmonised_orf_table
    removed_orf_ids = harmonise_orfs.out.removed_orf_ids

    emit:
    harmonised_orf_table
    removed_orf_ids  

}
