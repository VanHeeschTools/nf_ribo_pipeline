// Create MultiQC output report, including custom sections
process MULTIQC {
    label "multiqc"
    publishDir "${outdir}/multiqc/", mode: 'copy'

    input:
        path multiqc_files, stageAs: "?/*" // List, channel with paths to all files that should be in the report
        val multiqc_config                 // Path, multiqc_config 
        val outdir                         // Path, output directory

    output:
        path "*multiqc_report.html", emit: multiqc_report

    script:
        """
        multiqc . -c ${multiqc_config} -v
        """
}