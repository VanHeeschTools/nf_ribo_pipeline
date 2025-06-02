
process contaminants_check {

    tag "multi-sample-contaminants"
    label "samtools"
    publishDir "${outdir}/bowtie2", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(filtered_reads), path(sam_file)
    val keep_sam
    val outdir

    output:
    path "contaminant_counts_${meta.sample_id}_mqc.txt", emit: contaminant_samples
    path "passed_contaminant_counts_${meta.sample_id}_mqc.txt", emit: contaminant_samples_passed

    script:
    """
    sample_id="${meta.sample_id}"
    outfile="contaminant_counts_\${sample_id}_mqc.txt"
    outfile_passed="passed_contaminant_counts_\${sample_id}_mqc.txt"

    filtered_reads_n=\$(zcat "${filtered_reads}" | wc -l)
    filtered_reads_n=\$((filtered_reads_n / 4))

    read_counts=\$(samtools view -@ $task.cpus "${sam_file}" | awk '
        /rRNA/   {rRNA++}
        /tRNA/   {tRNA++}
        /snRNA/  {snRNA++}
        /snoRNA/ {snoRNA++}
        /mtDNA/  {mtDNA++}
        END {
          printf "%d\\t%d\\t%d\\t%d\\t%d", rRNA+0, tRNA+0, snRNA+0, snoRNA+0, mtDNA+0
        }
    ')
    echo -e "Sample\\tPassed\\trRNA\\ttRNA\\tsnRNA\\tsnoRNA\\tmtDNA" >> "\$outfile"
    echo -e "\$sample_id\\t\$filtered_reads_n\\t\$read_counts" >> "\$outfile"

    echo -e "Sample\\tPassed" >> "\$outfile_passed"
    echo -e "\$sample_id\\t\$filtered_reads_n" >> "\$outfile_passed"

    if [ "$keep_sam" = "false" ]; then
        rm -f "${sam_file}"
    fi
    """
}

process samtools {

    // Get mapping stats, sorted bam and .bai with SAMTOOLS

    tag "${meta.sample_id}"
    label "samtools"
    publishDir "${outdir}/star/", mode: 'copy'

    input: 
    tuple val(meta), path(bam) // Aligned BAMs
    val outdir                 // Output directory

    output:
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}*.Aligned.sortedByCoord.out.bam"), emit:sorted_bam
    path "${meta.sample_id}/${meta.sample_id}*.Aligned.sortedByCoord.out.bam", emit:bam_files
    path "${meta.sample_id}/${meta.sample_id}*" // Output all files to publishDir


    script:
    def new_bam = "${bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')}"
    def sample_id = meta.sample_id
    """
    mkdir -p ${sample_id}
    mkdir -p tmp/
    # Sort BAM
    samtools sort \
    -@ $task.cpus \
    -l 9 \
    -o "${sample_id}/${new_bam}" \
    -T "tmp/" \
    "${bam}"

    rm -r tmp/

    # Create mapping statistics with samtools
    # TODO: will this work for both local and end2end or will it overwrite itself?
    samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${sample_id}_stats.txt"

    # Index the bam with samtools
    samtools index -@ $task.cpus "${sample_id}/${new_bam}"
    """
}