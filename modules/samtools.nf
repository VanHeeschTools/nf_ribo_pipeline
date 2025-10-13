
// Create statistic output files of known contaminants in the input fastq file
process contaminants_check {

    tag "multi-sample-contaminants"
    label "samtools"
    publishDir "${outdir}/bowtie2/mqc_files", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(filtered_reads), val(bam_file)
    val keep_bam
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

    read_counts=\$(samtools view -@ $task.cpus "${bam_file}" | awk '
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

    if [ "$keep_bam" = false ]; then
        rm -f "${bam_file}"
    fi
    """
}

// Get mapping stats, sorted bam and .bai with SAMTOOLS
process samtools {

    tag "${sample_id}"
    label "samtools"
    publishDir "${outdir}/star/${sample_id}", mode: 'copy'

    input: 
    tuple val(sample_id), path(bam) // Aligned BAMs
    val outdir                      // Output directory

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}*.Aligned.sortedByCoord.out.bam"), emit:sorted_bam
    path "${sample_id}/${sample_id}*.Aligned.sortedByCoord.out.bam", emit:bam_files
    path "${sample_id}/${sample_id}*" // Output all files to publishDir

    script:
    def new_bam = "${bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')}"
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
    # TODO: will this work for both local and end2end or will it overwrite itself? Probably overwrite
    samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${new_bam}_stats.txt"

    # Index the bam with samtools
    samtools index -@ $task.cpus "${sample_id}/${new_bam}"
    """
}