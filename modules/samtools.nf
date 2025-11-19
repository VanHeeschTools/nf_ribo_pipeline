// Get mapping stats, sorted bam and .bai with SAMTOOLS
process samtools {

    tag "${sample_id}"
    label "samtools"
    publishDir "${outdir}/star/", mode: 'copy'

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