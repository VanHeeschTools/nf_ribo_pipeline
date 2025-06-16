process sample_psites {

  // Create per sample a BED file with sorted in-frame P-sites found in the sample

  tag "${meta.sample_id}"
  label "Ribo_Seq_R_scripts"
  publishDir "${outdir}/bedfiles", mode: 'copy'

  input:
  tuple val(meta), path(riboseqc_results)
  val outdir

  output:
  tuple val("${meta.sample_id}"), path("${meta.sample_id}_psites.sorted.bed"), emit: sample_psite_bed
  
  script:
    """
    psite_from_riboseqc.R ${riboseqc_results} 

    sort -T \$PWD -k1,1 -k2,2n "${meta.sample_id}_psites.bed" > "${meta.sample_id}_psites.sorted.bed"
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
  path "${orfcaller_gtf.baseName}.gtf_psites_p0.sorted.bed", emit: psite_bed

  script:
  //TODO: Find a better way to do this
  def extension = orfcaller_gtf.getName().tokenize('.')[-1]
  def fileType = extension == 'bed' ? 'bed' : extension == 'gtf' ? 'gtf' : 'Unknown'
  """
  inframe_psite_bed.py -i ${orfcaller_gtf} -a "no" -o \$PWD -t "${type}"

  sort -T \$PWD -k1,1 -k2,2n "${orfcaller_gtf.baseName}.${fileType}_psites_plus_partial.bed" > "${orfcaller_gtf.baseName}.gtf_psites_p0.sorted.bed"
  """

}


process merge_orfcaller_psites {
    label "merge_psites"

    input:
    path psite

    output:
    path "combined_psites.bed", emit: combined_psites

    script:
    """
    cat ${psite.join(' ')} | sort -T \$PWD --parallel=$task.cpus -k1,1 -k2,2n > combined_psites.bed
    """
}
