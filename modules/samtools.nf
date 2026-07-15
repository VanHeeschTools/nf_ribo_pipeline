// Get mapping stats, sorted bam and .bai with SAMTOOLS
process samtools {

    tag "${sample_id}"
    label "Ribo_Seq_tools"

    input: 
        tuple val(sample_id), path(bam) // Aligned BAMs

    output:
        tuple val(sample_id), path("${sample_id}/${sample_id}*.sortedByCoord.out.bam"), emit:sorted_bam
        path "${sample_id}/${sample_id}*.sortedByCoord.out.bam", emit:bam_files
        path "${sample_id}/${sample_id}*" // Output all files to publishDir

    when:
        task.ext.when == null || task.ext.when

    script:
        def new_bam = bam.name

        if (bam.name.endsWith(".Aligned.toTranscriptome.out.bam")) {
            new_bam = bam.name.replaceFirst('.Aligned.toTranscriptome.out.bam', '.Aligned.toTranscriptome.sortedByCoord.out.bam')
        } else if (bam.name.endsWith(".Aligned.out.bam")) {
            new_bam = bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')
        } else {
            error "Unexpected BAM filename: ${bam.name}"
        }

        //def new_bam = "${bam.name.replaceFirst('.Aligned.out.bam', '.Aligned.sortedByCoord.out.bam')}"
        """
        mkdir -p ${sample_id}
        mkdir -p tmp/
        # Sort BAM
        samtools sort \
        -@ $task.cpus \
        -o "${sample_id}/${new_bam}" \
        -T "tmp/" \
        "${bam}"

        rm -r tmp/

        # Create mapping statistics with samtools
        samtools stats -@ $task.cpus "${sample_id}/${new_bam}" > "${sample_id}/${new_bam}_stats.txt"

        # Index the bam with samtools
        samtools index -@ $task.cpus "${sample_id}/${new_bam}"
        """
}