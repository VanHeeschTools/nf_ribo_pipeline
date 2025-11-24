// Create index for STAR if a valid path is not given in parm file
process star_index {

    label "alignment"
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

    tag "${sample_id}"
    label "alignment"

    input: 
    tuple val(sample_id), path(reads)   // Trimmed RPF reads
    val outdir                          // Output directory
    val gtf                             // Transcriptome GTF file
    val star_index_path                 // STAR index

    output:
    path("${sample_id}/${sample_id}.*")
    tuple val(sample_id), path("${sample_id}/${sample_id}.local.Aligned.out.bam"), optional: true, emit: bams
    path "${sample_id}/${sample_id}.local.Log.final.out", emit: star_log_local

    script:
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

    tag "${sample_id}"
    label "alignment"
    publishDir "${outdir}/star/", mode: 'copy' , pattern: "${sample_id}/${sample_id}.end2end.Aligned.toTranscriptome.out.bam"


    input: 
    tuple val(sample_id), path(reads) // Trimmed RPF reads
    val outdir                        // Output directory
    val gtf                           // Transcriptome GTF file
    val star_index_path               // STAR index

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.end2end.Aligned.out.bam"), optional: true, emit: bams_end2end
    tuple val(sample_id), path("${sample_id}/${sample_id}.end2end.Aligned.toTranscriptome.out.bam"), optional: true, emit: bams_end2end_transcriptome
    path "${sample_id}/${sample_id}.end2end.Log.final.out", emit: star_log_end_to_end

    script:
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