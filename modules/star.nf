// Create index for STAR if a valid path is not given in parm file
process star_index {

    label "Ribo_Seq_tools"

    input: 
        val genome // Reference genome fasta file
        val gtf    // Transcriptome GTF file

    output:
        path "star_index", emit: star_index_path
        path "star_index/*"

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        STAR \
        --runMode genomeGenerate \
        --runThreadN $task.cpus \
        --sjdbGTFfile ${gtf} \
        --sjdbOverhang 29 \
        --genomeDir "star_index" \
        --genomeFastaFiles ${genome}
        """
}

// Aligns RPF reads to the reference genome
process star {
    tag "${sample_id}"
    label "Ribo_Seq_tools"
    label "alignment"
    
    input:
        tuple val(sample_id), path(reads) // Sample id and filtered reads
        val gtf                           // Reference gtf
        val star_index_path               // Path to star index directory
        val local                         // Boolean, true if running local mode
    
    output:
        path("${sample_id}/${sample_id}.*")
        tuple val(sample_id), path("${sample_id}/*.Aligned.out.bam"), emit: bam_file
        tuple val(sample_id), path("${sample_id}/*.Aligned.toTranscriptome.out.bam"), optional: true, emit: bam_file_transcriptome
        path "${sample_id}/*.Log.final.out", emit: star_log
    
    when:
        task.ext.when == null || task.ext.when
    
    script:
        def prefix = local ? "${sample_id}/${sample_id}.local." : "${sample_id}/${sample_id}.end2end."
        def extra_params = local ?
            "--outSAMattributes All" :
            "--alignEndsType EndToEnd --outSAMattributes MD NH --quantMode TranscriptomeSAM"
        """
        STAR \
            --genomeDir ${star_index_path} \\
            --sjdbGTFfile ${gtf} \\
            --readFilesIn ${reads} \\
            --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:ILLUMINA SM:${sample_id} \\
            --outFileNamePrefix ${prefix} \\
            --runThreadN ${task.cpus} \\
            --readFilesCommand zcat \\
            --outSAMtype BAM Unsorted \\
            --runDirPerm All_RWX \\
            --twopassMode Basic \\
            --outFilterMismatchNmax 2 \\
            --outFilterMultimapNmax 20 \\
            --outFilterType BySJout \\
            --alignSJoverhangMin 1000 \\
            --outTmpKeep None \\
            ${extra_params}
        """
}


// Run preseq
process preseq {
    tag "${sample_id}"
    label "Ribo_Seq_tools"

    input:
        tuple val(sample_id), path(sorted_bam_file) // Sample id and sorted local BAM file
        val c_curve // Bool true if c_curve should be run, otherwise lc_extrap will be run

    output:
        path "${out_file}", emit: preseq_txt

    script:
        out_file = c_curve ? "${sample_id}_c_curve.txt" : "${sample_id}_lc_extrap.txt"
        def subcommand = c_curve ? "c_curve" : "lc_extrap"
        def input_arguments = c_curve ?
            "-B -v -s 500000" :
            "-e 500000000 -s 1000000"

        """
        preseq ${subcommand} ${input_arguments} -o ${out_file} ${sorted_bam_file}
        """

}