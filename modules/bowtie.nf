process bowtie2_index {

    // Create index for bowtie2 alignment

    label "bowtie2"
    publishDir "${outdir}/bowtie2_index", mode: 'copy'

    input:
    path contaminants_fasta // FASTA with unwanted sequences
    val outdir              // Output directory

    output:
    path "bowtie2_index", emit: bowtie2_index_prefix
    path "bowtie2_index/bowtie2_index*.bt2"

    script:
    """
    mkdir -p "bowtie2_index"
    bowtie2-build \
    -f \
    --seed 24 \
    --threads $task.cpus \
    ${contaminants_fasta} \
    "bowtie2_index/bowtie2_index"
    """
}
 
process bowtie2 {

    // Use BOWTIE2 to align against unwanted sequences such
    // as rRNA, tRNA, snRNA, snoRNA, mtDNA and keep true
    // ribosome-protected fragments for mapping and ORF calling

    tag "${sample_id}"
    label "bowtie2"
    publishDir "${outdir}/bowtie2", mode: 'copy', pattern: "${sample_id}/${sample_id}_filtered.fastq.gz"


    input:
    val bowtie2_index_prefix      // Bowtie2 reference index
    tuple val(sample_id), path(reads)  // Trimmed reads
    val outdir                    // Output directory

    output:
    tuple val(sample_id), path(reads), path("${sample_id}/${sample_id}_filtered.fastq.gz"), path("${sample_id}/${sample_id}_contaminants.bam"), emit: bowtie_output_files
    tuple val(sample_id), path("${sample_id}/${sample_id}_filtered.fastq.gz"), emit: filtered_reads

    script:
    """
    mkdir -p "${sample_id}"
    bowtie2 \
    --seedlen=25 \
    --threads $task.cpus \
    --time \
    --un-gz "${sample_id}/${sample_id}_filtered.fastq.gz" \
    -x ${bowtie2_index_prefix} \
    -U ${reads} |samtools view -@ $task.cpus -bS - > "${sample_id}/${sample_id}_contaminants.bam"
    """
}


// Create statistic output files of known contaminants in the input fastq file
process contaminants_check {

    tag "multi-sample-contaminants"
    label "samtools"
    publishDir "${outdir}/bowtie2/mqc_files", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(filtered_reads), val(bam_file)
    val keep_bam
    val outdir

    output:
    path "contaminant_counts_${sample_id}_mqc.txt", emit: contaminant_samples
    path "passed_contaminant_counts_${sample_id}_mqc.txt", emit: contaminant_samples_passed

    script:
    """
    sample_id="${sample_id}"
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