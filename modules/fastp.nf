// Remove low-quality and short reads from fastq file using fastp
// Currently not being used in favor of trimgalore

process read_trimming {
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

// Remove low-quality and short reads from fastq file using trimgalore
process trimgalore{

    tag "${meta.sample_id}"
    label "read_trimming"

    input:
    tuple val(meta),val(reads)    // Input FASTQ reads
    val outdir                    // Output directory

    output:
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}_trimmed.fq.gz"), emit: reads
    path "${meta.sample_id}/${meta.sample_id}_trimming_report.txt", emit: trimgalore_trimming_report
    path "total_reads_${meta.sample_id}_mqc.txt", emit: total_reads
    path "removed_reads_${meta.sample_id}_mqc.txt", emit: removed_reads

    script:
    """
    # Obtain total number of reads at the start and write to multiqc supported file format
    outfile_total="total_reads_${meta.sample_id}_mqc.txt"
    outfile_removed="removed_reads_${meta.sample_id}_mqc.txt"
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

    # Rename trimming report so it mathces the sample_id
    mv ${meta.sample_id}/*_trimming_report.txt ${meta.sample_id}/${meta.sample_id}_trimming_report.txt

    # Obtain amount of removed reads
    num=\$(grep "Sequences removed" ${meta.sample_id}/${meta.sample_id}_trimming_report.txt | awk -F':' '{print \$2}' | awk '{print \$1}')

    echo -e "Sample\\tTotal input reads\\tRemoved" >> "\$outfile_removed"
    echo -e "${meta.sample_id}\\t\${total_reads_n}\\t\${num}" >> "\$outfile_removed"
    """

}
