process bowtie2_index {

    // Create index for bowtie2 alignment

    label "bowtie2_index"
    publishDir "${outdir}/bowtie2_index", mode: 'copy'

    input:
    path contaminants_fasta // FASTA with unwanted sequences
    val outdir             // Output directory

    output:
    val "bowtie2_index" , emit: bowtie2_index_prefix
    path "bowtie2_index*.bt2"

    script:
    """
    bowtie2-build \
    -f \
    --seed 24 \
    --threads $task.cpus \
    ${contaminants_fasta} \
    bowtie2_index

    """
}
 
process bowtie2 {

    // Use BOWTIE2 to align against unwanted sequences such
    // as rRNA, tRNA, snRNA, snoRNA, mtDNA and keep true
    // ribosome-protected fragments for mapping and ORF calling

    tag "${meta.sample_id}"
    label "bowtie2"
    publishDir "${outdir}/bowtie2", mode: 'copy'

    input:
    val bowtie2_index_prefix      // Bowtie2 reference index
    tuple val(meta), path(reads)  // Trimmed fastp reads
    val outdir                    // Output directory

    output:
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}_filtered.fastq.gz"), emit: filtered_reads
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}_contaminants.sam"), emit: sam_file

    script:
    def sample_id = meta.sample_id
    """
    mkdir "${sample_id}"
    bowtie2 \
        --seedlen=25 \
        --threads $task.cpus \
        --time \
        --un-gz "${sample_id}/${sample_id}_filtered.fastq.gz" \
        -x ${bowtie2_index_prefix} \
        -U ${reads} \
        -S "${sample_id}/${sample_id}_contaminants.sam"
    """
}
