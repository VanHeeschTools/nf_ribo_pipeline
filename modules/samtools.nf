
process contaminants_check {

    // Use SAMTOOLS to output the reads per contaminant for QC

    tag "${meta.sample_id}"
    label "contaminants"
    publishDir "${outdir}/bowtie2", mode: 'copy'

    input:
    tuple val(meta), path(reads)            // Trimmed fastp reads
    tuple val(meta_2), path(filtered_reads)   // RPF reads
    tuple val(meta_3), path(sam_file)         // Filtered contaminant reads
    val bowtie2_index_prefix                // Bowtie2 reference index
    val outdir                              // Output directory
    val keep_sam                            // Boolean, keep big SAM file for debugging

    output:

    path "${meta.sample_id}/${meta.sample_id}_*.txt", emit: contaminants_bowtie2

    script:
    def sample_id = meta.sample_id
    """
    #!/bin/bash

    mkdir ${sample_id}

    # Create contaminant QC file
    contaminants_type=\$(basename "${bowtie2_index_prefix}")
    contaminants_file="${sample_id}/${sample_id}_\${contaminants_type}.txt"

    touch "\${contaminants_file}"

    # Get total number of reads
    tot_reads=\$(zcat "${reads}" | wc -l)
    tot_reads=\$((tot_reads / 4))

    echo -e "RiboseQC run for ${sample_id} on \$(date) \n" >> "\${contaminants_file}"

    # Print headers to file
    printf '\t%s\t%s\t%s\n' "READ_TYPE" "READS" "PERCENTAGE" >> "\${contaminants_file}"

    # Print total no. of reads
    printf '%s\t%s\t%s\t%s\n' "${sample_id}" "Total" "\${tot_reads}" "100" >> "\${contaminants_file}"

    # For each contaminant type, print absolute and relative number of reads
    for contaminant_type in tRNA snRNA snoRNA mtDNA rRNA; do
        contaminant_reads=\$(samtools view -@ $task.cpus "${sam_file}" | grep -c "\${contaminant_type}")
        contaminant_percentage=\$(awk -v n=\${tot_reads} -v r=\${contaminant_reads} 'BEGIN{printf "%.2f", r/n*100}')
        printf '%s\t%s\t%s\t%s\n' "${sample_id}" "\${contaminant_type}" "\${contaminant_reads}" "\${contaminant_percentage}" >> "\${contaminants_file}"
    done

    # Count reads that passed filtering
    filtered_reads_n=\$(zcat "${filtered_reads}" | wc -l)
    filtered_reads_n=\$((filtered_reads_n / 4))
    filtered_percentage=\$(awk -v n=\${tot_reads} -v r=\${filtered_reads_n} 'BEGIN{printf "%.2f", r/n*100}')
    printf '%s\t%s\t%s\t%s\n\n' "${sample_id}" "Passed" "\${filtered_reads_n}" "\${filtered_percentage}" >> "\${contaminants_file}"

    # Remove the SAM file if keep_sam is false
    if [ "${keep_sam}" = "false" ]; then
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
    val outdir                      // Output directory

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
    samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${sample_id}_stats.txt"

    # Index the bam with samtools
    samtools index -@ $task.cpus "${sample_id}/${new_bam}"
    """
}
