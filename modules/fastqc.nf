// Generate QC files using fastqc, shown in the MultiQC report
process fastqc {

    tag "${sample_id}"
    label "fastqc"
    publishDir "${outdir}/fastqc", mode: 'copy'

    input:
        tuple val(sample_id), path(reads) // Tuple, meta info plus trimmed FASTQ reads
        val outdir                        // Path, output directory

    output:
        path "${sample_id}/${sample_id}_filtered_fastqc.html", emit: fastqc_html // Output QC summary 
        path "${sample_id}/${sample_id}_filtered_fastqc.zip",  emit: fastqc_zip  // QC files

    script:
        """
        # Create temp directory to run fastqc
        mkdir -p tmp
        mkdir -p ${sample_id}

        # Run fastqc
        fastqc \
        ${reads} \
        --threads $task.cpus \
        -d "tmp" \
        --outdir "${sample_id}" 

        # Remove fastqc temp direcory
        rm -r tmp
        """
}