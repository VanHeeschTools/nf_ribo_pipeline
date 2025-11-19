// Create RiboTIE template to be used in the next steps
process create_template {
    label "create_ribotie_template"
    publishDir "${outdir}/ribotie", mode: 'copy'

    input:
    val ribotie_bams // list of tuples of sample id and bam path
    path gtf
    path fasta
    val outdir

    output:
    path "ribotie_template.yml", emit: template

    //TODO: fix starting from this point
    script:
    def echo_lines = ribotie_bams.collect { pair ->
        def sample = pair[0]
        def bam = pair[1]
        // echo one line per sample
        "echo \"  ${sample} : ${bam}\" >> ribotie_template.yml"
    }.join('\n')

    """
    echo "gtf_path : ${gtf}" > ribotie_template.yml
    echo "fa_path : ${fasta}" >> ribotie_template.yml
    echo "ribo_paths :" >> ribotie_template.yml
    ${echo_lines}
    """
}

// Create h5 database based on genomic features from the gtf 
process parse_genomic_features {
    label "ribotie"
    publishDir "${outdir}/ribotie", mode: 'copy'

    input:
    path gtf
    path fasta
    val outdir

    output:
    path "genomic_features_db.h5", emit: h5_path

    script:
    """
    tis_transformer \
    --gtf_path ${gtf} \
    --fa_path ${fasta} \
    --h5_path genomic_features_db.h5 \
    --data  \
    --no_backup \
    --num_workers $task.cpus 
    """
}

// Create h5 database for every sample
process parse_samples {
    label "ribotie"
    publishDir "${outdir}/ribotie", mode: 'copy'

    input:
    path genomic_h5_db
    path ribotie_template
    tuple val(sample_id), path(sample_bam)
    path gtf
    path fasta
    val outdir

    output:
    tuple val(sample_id), path("genomic_features_db_${sample_id}.h5"), emit: h5_path

    script:
    """
    ribotie ${ribotie_template} \
    --h5_path ${genomic_h5_db} \
    --data \
    --samples ${sample_id} \
    --parallel \
    --no_backup \
    --num_workers $task.cpus \
    --accelerator cpu
    """
}

// Run RiboTIE for all samples individually
process ribotie_predict_samples {
    label "ribotie"
    publishDir "${outdir}/ribotie", mode: 'copy'

    input:
    tuple val(sample_id), path(sample_h5)
    path genomic_h5_db
    path ribotie_template
    path gtf
    path fasta
    val outdir

    output:
    path "genomic_features_db_${sample_id}.csv", emit: ribotie_orf_csv
    path "genomic_features_db_${sample_id}.gtf", emit: ribotie_orf_gtf

    script:
    """
    ribotie ${ribotie_template} \
    --h5_path ${genomic_h5_db} \
    --samples ${sample_id} \
    --parallel \
    --num_workers $task.cpus \
    --no_backup \
    --return_ORF_coords

    mv multiqc multiqc_${sample_id} # Rename to make unique
    """
}

// Merge and filter the RiboTIE output files
process merge_ribotie_output{
    publishDir "${outdir}/merged_ribotie", mode: 'copy', pattern: "*.csv"
    label "ribotie"

    input:
    val ribotie_csv_files
    path genomic_h5_db
    val ribotie_min_samples
    val outdir

    output:
    path "RiboTIE_merged.gtf", emit: ribotie_merged_gtf
    path "RiboTIE_merged.csv"
    path "RiboTIE_duplicate_filtered_merged.csv"
    path "RiboTIE_unfiltered_merged.csv"

    script:
    """    
    merge_ribotie.py ${genomic_h5_db} "${ribotie_csv_files.join(',')}" ${ribotie_min_samples}
    """
}

process ribotie_add_stop{
    publishDir "${outdir}/merged_ribotie", mode: 'copy'
    label "Ribo_Seq_R_scripts"

    input:
    path merged_ribotie
    path gtf
    val outdir

    output:
    path "RiboTIE.gtf", emit: ribotie_gtf

    script:
    """    
    fix_ribotie_stop.R ${merged_ribotie} ${gtf}
    """
}

