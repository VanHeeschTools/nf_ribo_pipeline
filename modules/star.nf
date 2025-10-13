// Create index for STAR if a valid path is not given in parm file
process star_index {

    label "alingment"
    publishDir "${outdir}/star_index/", mode: 'copy'

    input: 
    val genome // Reference genome fasta file
    val gtf    // Transcriptome GTF file
    val outdir // Output directory

    output:
    path "star_index", emit: star_index_path
    path "star_index/*"

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

// Aligns RPF reads to the reference genome to create input for RiboseQC
process star_local{

    tag "${meta.sample_id}"
    label "alignment"

    input: 
    tuple val(meta), path(reads)   // Trimmed RPF reads
    val outdir                     // Output directory
    val gtf                        // Transcriptome GTF file
    val star_index_path            // STAR index

    output:
    path("${meta.sample_id}/${meta.sample_id}.*")
    tuple val(meta.sample_id), path("${meta.sample_id}/${meta.sample_id}.local.Aligned.out.bam"), optional: true, emit: bams
    path "${meta.sample_id}/${meta.sample_id}.local.Log.final.out", emit: star_log_local

    script:
    def sample_id = meta.sample_id

    """
    # ORFquant BAM
    STAR \
    --genomeDir ${star_index_path} \
    --sjdbGTFfile ${gtf} \
    --readFilesIn ${reads} \
    --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id} \
    --outFileNamePrefix "${sample_id}/${sample_id}.local." \
    --runThreadN $task.cpus \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --runDirPerm All_RWX \
    --twopassMode Basic \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 20 \
    --outSAMattributes All \
    --outFilterType BySJout \
    --alignSJoverhangMin 1000 \
    --outTmpKeep None
    """

}

// Aligns RPF reads to the reference genome to create PRICE input
process star_end_to_end {

    tag "${meta.sample_id}"
    label "alignment"

    input: 
    tuple val(meta), path(reads) // Trimmed RPF reads
    val outdir                   // Output directory
    val gtf                      // Transcriptome GTF file
    val star_index_path          // STAR index

    output:
    tuple val(meta.sample_id), path("${meta.sample_id}/${meta.sample_id}.end2end.Aligned.out.bam"), optional: true, emit: bams_end2end
    tuple val(meta.sample_id), path("${meta.sample_id}/${meta.sample_id}.end2end.Aligned.toTranscriptome.out.bam"), optional: true, emit: bams_end2end_transcriptome
    path "${meta.sample_id}/${meta.sample_id}.end2end.Log.final.out", emit: star_log_end_to_end

    script:
    def sample_id = meta.sample_id

    """
    STAR \
    --genomeDir ${star_index_path} \
    --sjdbGTFfile ${gtf} \
    --readFilesIn ${reads} \
    --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id} \
    --outFileNamePrefix "${sample_id}/${sample_id}.end2end." \
    --runThreadN $task.cpus \
    --quantMode TranscriptomeSAM \
    --readFilesCommand zcat \
    --outSAMtype BAM Unsorted \
    --runDirPerm All_RWX \
    --twopassMode Basic \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 20 \
    --outSAMattributes MD NH \
    --outFilterType BySJout \
    --alignSJoverhangMin 1000 \
    --alignEndsType EndToEnd \
    --outTmpKeep None
    """
}