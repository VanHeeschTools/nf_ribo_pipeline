 // Create per sample a BED file with sorted in-frame P-sites found in the sample
process sample_psites {

    tag "${sample_id}"
    label "Ribo_Seq_R_scripts"

    input:
      tuple val(sample_id), path(riboseqc_results)

    output:
      tuple val("${sample_id}"), path("${sample_id}_psites.sorted.bed"), emit: sample_psite_bed

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      psite_from_riboseqc.R ${riboseqc_results} 
      sort -T \$PWD -k1,1 -k2,2n "${sample_id}_psites.bed" > "${sample_id}_psites.sorted.bed"
      """
}

// Create reference in-frame P-sites file from GTF
process reference_psites {

    label "Ribo_Seq_R_scripts"

    input:
      path orfcaller_gtf
      val type
      path reference_protein_fa
      path package_install_loc

    output:
      path "${orfcaller_gtf.baseName}_p0_reference_sorted.bed", emit: reference_psite_bed
      path "${orfcaller_gtf.baseName}_correct_cds.rds",         emit: reference_cds_rds

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      create_p0_bed.R ${orfcaller_gtf} ${type} ${reference_protein_fa} ${package_install_loc}
      sort -T \$PWD -k1,1 -k2,2n "${orfcaller_gtf.baseName}_p0.bed" > "${orfcaller_gtf.baseName}_p0_reference_sorted.bed"
      """

}

// Create reference in-frame P-sites file from GTF
process orfcaller_psites {

    label "Ribo_Seq_R_scripts"

    input:
      path orfcaller_gtf 
      val type
      path reference_protein_fa
      path package_install_loc

    output:
      tuple path(orfcaller_gtf), path("${orfcaller_gtf.baseName}_p0_orf_sorted.bed"), emit: orfcaller_psite_bed
      path "${orfcaller_gtf.baseName}_p0_orf_sorted.bed", emit: orf_psite_bed

    when:
        task.ext.when == null || task.ext.when


    script:
      """
      create_p0_bed.R ${orfcaller_gtf} ${type} ${reference_protein_fa} ${package_install_loc}
      sort -T \$PWD -k1,1 -k2,2n "${orfcaller_gtf.baseName}_p0.bed" > "${orfcaller_gtf.baseName}_p0_orf_sorted.bed"
      """

}

// Combine the p_site bed files from all ORFcallers
process merge_orfcaller_psites {

    label "merge_psites"

    input:
      path orfcaller_psites

    output:
      path "combined_psites_unfiltered.bed", emit: combined_psites

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      cat ${orfcaller_psites.join(' ')} | sort -T \$PWD --parallel=$task.cpus -k1,1 -k2,2n > combined_psites_unfiltered.bed
      """
}

process orf_ref_p0_intersect {

    label "intersect_psites"

    input:
      tuple path(orfcaller_gtf), path(orf_psite_bed)
      path ref_psite_bed

    output:
      tuple path(orfcaller_gtf), path("${orfcaller_gtf.baseName}_ref_intersect.bed"), emit: orf_ref_intersect

    when:
        task.ext.when == null || task.ext.when

    script:
    """
      bedtools intersect \
      -a ${ref_psite_bed} \
      -b ${orf_psite_bed} \
      -wa \
      -wb \
      -header \
      -f 1.00 \
      -s \
      -sorted > "${orfcaller_gtf.baseName}_ref_intersect.bed"
      """
}