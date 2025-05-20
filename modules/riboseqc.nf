process create_annotation {

    label "orfquant_annotation"
    publishDir "${outdir}/orfquant_annotation", mode: 'copy'

    input:
    val gtf
    val twobit
    val package_install_loc
    val orfquant_prefix

    output:
    val "${orfquant_prefix}_Rannot", emit: orfquant_annotation

    script:
    """
    create_riboseq_annotation.R \
    ${twobit} \
    ${gtf} \
    ${package_install_loc} \
    ${orfquant_prefix}
    """

}

process riboseqc {

    tag "${meta.sample_id}"
    label "riboseqc"
    publishDir "${outdir}/riboseqc", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    val outdir
    val orfquant_annotation
    val pandoc_dir
    val orfquant_annot_package
    val package_install_loc

    output:
    tuple val(meta),path("${meta.sample_id}/${meta.sample_id}_for_ORFquant"), emit: orfquant_psites
    path("${meta.sample_id}/${meta.sample_id}*")

    script:
    """
    run_riboseqc.R \
        ${bam} \
        ${meta.sample_id}/${meta.sample_id} \
        ${orfquant_annotation} \
        ${pandoc_dir} \
        ${orfquant_annot_package} \
        ${package_install_loc}
    """
}
