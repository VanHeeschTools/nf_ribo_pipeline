// Generate QC files using fastqc, shown in the MultiQC report
process fastqc {

    tag "${sample_id}"
    label "Ribo_Seq_tools"

    input:
        tuple val(sample_id), path(reads) // Tuple, meta info plus trimmed FASTQ reads

    output:
        path "${sample_id}_filtered_fastqc.html", emit: fastqc_html // Output QC summary 
        path "${sample_id}_filtered_fastqc.zip",  emit: fastqc_zip  // QC files

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        # Create temp directory to run fastqc
        mkdir -p tmp
        mkdir -p ${sample_id}

        # Run fastqc
        fastqc \
        ${reads} \
        --threads $task.cpus \
        --dir "tmp" \
        --outdir "." 

        # Remove fastqc temp direcory
        rm -r tmp
        """
}