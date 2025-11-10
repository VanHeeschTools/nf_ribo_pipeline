// Generate QC files using fastqc, shown in the MultiQC report
process fastqc {

    tag "${meta.sample_id}"
    label "fastqc"
    publishDir "${outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads) // Tuple, meata info plus trimmed FASTQ reads
    val outdir                   // Path, output directory

    output:
    path "${meta.sample_id}/${meta.sample_id}_filtered_fastqc.html" // Output QC summary 
    path "${meta.sample_id}/${meta.sample_id}_filtered_fastqc.zip", emit: fastqc_zip  // QC files

    script:
    def sample_id = meta.sample_id
    """
    mkdir -p tmp
    mkdir -p ${sample_id}
    fastqc \
    ${reads} \
    --threads $task.cpus \
    -d "tmp" \
    --outdir "${sample_id}" 
    rm -r tmp
    """

}

// Generate size_distibution between 20-40 nucleotides 
process size_distribution {
    label "Ribo_Seq_R_scripts"

    input:
    tuple val(meta), path(reads) // Tuple, meata info plus trimmed FASTQ reads

    output:
    path "${meta.sample_id}_size_distribution_mqc.txt", emit: size_distribution

    script:
    """
    size_distribution.R ${meta.sample_id} ${reads}
    """

}
