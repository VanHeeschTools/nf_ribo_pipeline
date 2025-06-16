process MULTIQC {
    label "multiqc"
    publishDir "${outdir}/multiqc/", mode: 'copy'

    input:
    path multiqc_files, stageAs: "?/*" // Channel with paths to all files that should be in the report
    val multiqc_config                 // Path of the multiqc_config 
    val outdir                         // Path of output directory

    output:
    path "*multiqc_report.html", emit: multiqc_report

    script:
    """
    multiqc . -c ${multiqc_config} -v
    """
}