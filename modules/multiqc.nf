// Create MultiQC output report, including custom sections
process MULTIQC {
    label "multiqc"

    input:
        path multiqc_files, stageAs: "?/*" // List, channel with paths to all files that should be in the report
        val multiqc_config                 // Path, multiqc_config 

    output:
        path "*multiqc_report.html", emit: multiqc_report

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        multiqc . -c ${multiqc_config} -v
        """
}