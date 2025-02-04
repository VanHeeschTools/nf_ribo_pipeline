process ref_psites {

    // Create reference in-frame P-sites file from GTF

    label "ref_psites"
    publishDir "${outdir}/annotation", mode: 'copy'

    input:
    val gtf
    val outdir

    output:
    path(${ref_base}.gtf_psites_p0.sorted.bed), emit: ref_psite_bed

    script:
    // return basename for second part of the script
    def ref_base = gtf.baseName.replaceAll(/\.gtf$/, '')
    """
    python3 calculate_psite_coords_final.py \
    -i ${gtf} \
    -a "no" \
    -o `pwd` \
    -t "ORF_id"

    sort -k1,1 -k2,2n "${ref_base}.gtf_psites_plus_partial.bed" > "${ref_base}.gtf_psites_p0.sorted.bed"
    """

}

process sample_psites {

    // Create per sample a BED file with sorted in-frame P-sites found in the sample

    tag "${meta.sample_id}"
    label "sample_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

  input:
  tuple val(meta), val(riboseqc_results)
  val package_install_loc
  val outdir

  output:
  tuple val(meta), path("${meta.sample_id}_psites.sorted.bed"), emit: sample_psite_bed
  
  script:
    """
    Rscript psites_bed_from_riboseqc.R \
    ${riboseqc_results} 
    ${package_install_loc}
    sort -k1,1 -k2,2n "${meta.sample_id}_psites.bed" > "${meta.sample_id}_psites.sorted.bed"
    """
}
