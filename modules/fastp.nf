process read_trimming {

    // Remove low-quality and short reads from fastq file

    tag "${meta.sample_id}"
    label "read_trimming"

    input:
    tuple val(meta),val(reads)    // Input FASTQ reads
    val outdir                    // Output directory

    output:
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}.fastp.fastq.gz") , emit: reads
    path("${meta.sample_id}/${meta.sample_id}.json")
    path("${meta.sample_id}/${meta.sample_id}.html")
    path("${meta.sample_id}/${meta.sample_id}.log")

    script:
    def sample_id = meta.sample_id
    //TODO: compare with trimgalore output
    """
    mkdir "${sample_id}"
    fastp \
    --thread $task.cpus \
    --in1 "${reads}" \
    --out1 "${sample_id}/${sample_id}.fastp.fastq.gz" \
    --json ${sample_id}/${sample_id}.json \
    --html ${sample_id}/${sample_id}.html \
    --length_required 25 \
    2> >(tee ${sample_id}/${sample_id}.log >&2)
    """
}

process trimgalore{

    tag "${meta.sample_id}"
    label "read_trimming"

    input:
    tuple val(meta),val(reads)    // Input FASTQ reads
    val outdir                    // Output directory
 
    output:
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}_trimmed.fq.gz"), emit: reads
    path "total_reads_${meta.sample_id}_mqc.txt", emit: total_reads_samples

    script:
    """
    # Obtain total number of reads at the start and write to multiqc supported file format
    outfile_total="total_reads_${meta.sample_id}_mqc.txt"
    total_reads_n=\$(zcat "${reads}" | wc -l)
    total_reads_n=\$((total_reads_n / 4))
    echo -e "Sample\\tTotal" >> "\$outfile_total"
    echo -e "${meta.sample_id}\\t\$total_reads_n" >> "\$outfile_total"

    mkdir "${meta.sample_id}"
    
    trim_galore \
    "${reads}" \
    --cores $task.cpus \
    --gzip \
    --length 25 \
    --trim-n \
    --output_dir "${meta.sample_id}/" \
    --basename "${meta.sample_id}"
    """

}
