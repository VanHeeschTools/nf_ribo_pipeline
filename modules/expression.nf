process intersect_psites {

    // Intersect a reference BED with p-site positions with p-sites from a sample

    tag "${sample_id}"
    label "intersect_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_psite_bed)
    path ref_psite_bed
    val orfcaller
    val outdir

    output:
    tuple val(sample_id), path("${sample_id}_${orfcaller}_intersect.bed"), emit: sample_intersect_bed

    script:
    """
    echo ${sample_id}
    echo ${sample_psite_bed}
    echo ${ref_psite_bed}
    echo ${outdir}


    bedtools intersect \
      -a ${sample_psite_bed} \
      -b ${ref_psite_bed} \
      -wa \
      -wb \
      -header \
      -f 1.00 \
      -s \
      -sorted > "${sample_id}_${orfcaller}_intersect.bed"
      """
}


process ppm_matrix {

    // Create a matrix object for raw P-sites and PPM for each ORF and
    // each sample included in the cohort

    label "calculate_matrix"
    publishDir "${outdir}/orf_expression", mode: 'copy'

    input:
    path ref_psite_bed
    path sample_intersect_bed
    val orfcaller_name
    val outdir

    output:
    path("${orfcaller_name}_psites_permillion.csv"), emit: ppm_matrix
    path("${orfcaller_name}_psites.csv"), emit: psite_matrix

    script:
    """
    psite_matrix.R \
    "${ref_psite_bed}" \
    "${sample_intersect_bed}" \
    "${orfcaller_name}"
    """
}
