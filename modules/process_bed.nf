 // Create per sample a BED file with sorted in-frame P-sites found in the sample
process sample_psites {

    tag "${sample_id}"
    label "Ribo_Seq_R"

    input:
      tuple val(sample_id), path(riboseqc_results)

    output:
      tuple val("${sample_id}"), path("${sample_id}_psites.sorted.bed.gz"), emit: sample_psite_bed
      path "${sample_id}_psites.sorted.bed.gz"

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      psite_from_riboseqc.R ${riboseqc_results} 
      sort -T \$PWD -k1,1 -k2,2n "${sample_id}_psites.bed" > "${sample_id}_psites.sorted.bed"
      gzip -f "${sample_id}_psites.sorted.bed"
      """
}

// Create reference in-frame P-sites file from GTF
process orfcaller_psites {

    label "Ribo_Seq_R"

    input:
      path orfcaller_gtf 
      val type
      path package_install_loc

    output:
      tuple path(orfcaller_gtf), path("${orfcaller_gtf.baseName}_p0_orf_sorted.bed.gz"), emit: orfcaller_psite_bed
      path "${orfcaller_gtf.baseName}_p0_orf_sorted.bed.gz", emit: orf_psite_bed

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      process_ORF_and_Ref_gtf.R ${orfcaller_gtf} ${type} "null" ${package_install_loc}
      sort -T \$PWD -k1,1 -k2,2n "${orfcaller_gtf.baseName}_p0.bed" | gzip > "${orfcaller_gtf.baseName}_p0_orf_sorted.bed.gz"
      """
}

// Combine the p_site bed files from all ORFcallers
process merge_orfcaller_psites {

    label "Ribo_Seq_tools"

    input:
      path orfcaller_psites

    output:
      path "combined_psites_unfiltered.bed", emit: combined_psites

    when:
        task.ext.when == null || task.ext.when

    script:
      """
      zcat ${orfcaller_psites.join(' ')} | sort -T \$PWD --parallel=$task.cpus -k1,1 -k2,2n > combined_psites_unfiltered.bed
      gzip -kf combined_psites_unfiltered.bed
      """
}

