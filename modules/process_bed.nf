process sample_psites {

  // Create per sample a BED file with sorted in-frame P-sites found in the sample

  tag "${meta.sample_id}"
  label "sample_psites"
  publishDir "${outdir}/bedfiles", mode: 'copy'

  input:
  tuple val(meta), path(riboseqc_results)
  path package_install_loc
  val outdir

  output:
  tuple val("${meta.sample_id}"), path("${meta.sample_id}_psites.sorted.bed"), emit: sample_psite_bed
  
  script:
    """
    psite_from_riboseqc.R \
    ${riboseqc_results} \
    ${package_install_loc}

    sort -T \$PWD -k1,1 -k2,2n "${meta.sample_id}_psites.bed" > "${meta.sample_id}_psites.sorted.bed"
    """
}


process ref_psites {

   // Create reference in-frame P-sites file from GTF

  label "ref_psites"
  publishDir "${outdir}/annotation", mode: 'copy'

  input:
  path gtf
  val outdir

  output:
  path "${gtf.baseName}.gtf_psites_p0.sorted.bed", emit: ref_psite_bed

  script:
  """
  inframe_psite_bed.py -i ${gtf} -a "no" -o \$PWD -t "ORF_id"

  sort -T \$PWD -k1,1 -k2,2n "${gtf.baseName}.bed_psites_plus_partial.bed" > "${gtf.baseName}.gtf_psites_p0.sorted.bed"
  """

}