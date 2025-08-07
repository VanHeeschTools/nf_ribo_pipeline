process sample_psites {

  // Create per sample a BED file with sorted in-frame P-sites found in the sample

  tag "${sample_id}"
  label "Ribo_Seq_R_scripts"
  publishDir "${outdir}/bedfiles", mode: 'copy'

  input:
  tuple val(sample_id), path(riboseqc_results)
  val outdir

  output:
  tuple val("${sample_id}"), path("${sample_id}_psites.sorted.bed"), emit: sample_psite_bed

  script:
  """
  psite_from_riboseqc.R ${riboseqc_results} 
  sort -T \$PWD -k1,1 -k2,2n "${sample_id}_psites.bed" > "${sample_id}_psites.sorted.bed"
  """
}

process orfcaller_psites {

  // Create reference in-frame P-sites file from GTF

  label "obtain_psites"
  publishDir "${outdir}/annotation", mode: 'copy'

  input:
  path orfcaller_gtf
  val type
  val outdir

  output:
  path "${orfcaller_gtf.baseName}_gtf_psites_p0_sorted.bed", emit: psite_bed

  script:
  """
  inframe_psite_bed.py -i ${orfcaller_gtf} -a "no" -o \$PWD -t "${type}"
  sort -T \$PWD -k1,1 -k2,2n "${orfcaller_gtf.baseName}_psites_plus_partial.bed" > "${orfcaller_gtf.baseName}_gtf_psites_p0_sorted.bed"
  """

}

process merge_orfcaller_psites {

  // Combine the p_site bed files from all ORFcallers

  label "merge_psites"

  input:
  path orfcaller_psites

  output:
  path "combined_psites.bed", emit: combined_psites

  script:
  """
  cat ${orfcaller_psites.join(' ')} | sort -T \$PWD --parallel=$task.cpus -k1,1 -k2,2n |uniq > combined_psites.bed
  """
}
