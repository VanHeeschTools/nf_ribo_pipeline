include { orfcaller_psites; sample_psites; orfcaller_psites as ref_psites; merge_orfcaller_psites } from "../modules/process_bed.nf"

workflow PSITE {

    take:
    orfcaller_gtf
    reference_gtf
    riboseqc_results
    package_install_loc
    outdir

    main:
    
    // Create reference in-frame bed file for the ORF caller
    orfcaller_psites(orfcaller_gtf,
               "ORF_id",
               outdir
    )
   
    // Merge the bed files of all ORFcallers
    merge_orfcaller_psites(orfcaller_psites.out.collect()
    )

    // Create reference in-frame bed file for the reference gtf
    ref_psites(reference_gtf,
                "transcript_id",
                outdir)

    orfcaller_psites = merge_orfcaller_psites.out.combined_psites
    ref_psites = ref_psites.out.psite_bed

    emit:
    orfcaller_psites
    ref_psites

}
