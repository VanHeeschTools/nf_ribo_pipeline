include { annotate_orfs; harmonise_orfs } from "../modules/annotation.nf"

workflow ANNOTATION {

    take:
    orfcaller_output
    orfcaller_psites
    ref_psites
    reference_gtf
    package_install_loc
    orfquant_annot_package
    outdir
   
    main:
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
    // TODO: Maybe find a better way to do this, this seems risky
    harmonise_input = annotate_orfs.out.basic_orf_table
    .collect()
    .map { files ->
        def sorted = files.sort { it.name }
        tuple( sorted[0], sorted[1] )
    }

    harmonise_orfs(harmonise_input,
                   outdir)

    harmonised_orf_table = harmonise_orfs.out.harmonised_orf_table
    removed_orf_ids = harmonise_orfs.out.removed_orf_ids


    emit:
    harmonised_orf_table
    removed_orf_ids 

}
