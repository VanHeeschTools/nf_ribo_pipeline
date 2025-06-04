process fastqc {

    // Generate QC files

    tag "${meta.sample_id}"
    label "fastqc"
    publishDir "${outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads) // Trimmed FASTQ read from fastp
    val outdir                   // Output directory

    output:
    path "*/*.html" // Output QC summary 
    path "*/*.zip", emit: fastqc_zip  // QC files

    script:
    //TODO: put in config file
    def sample_id = meta.sample_id
    """
    mkdir ${sample_id}
    fastqc \
        ${reads} \
        --threads $task.cpus \
        --outdir "${sample_id}"
    """

}
