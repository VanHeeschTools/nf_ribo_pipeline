// Run RiboseQC index file creation
process riboseqc_index {

    label "Ribo_Seq_R"

    input:
        path reference_twobit
        path reference_gtf
        path reference_fasta

    output:
        path "bsgenome_install", emit: bsgenome_install_dir
        path "${reference_gtf}_Rannot", emit: rannot_file

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        mkdir -p bsgenome_install
        mkdir -p temp
        export TMPDIR=\$PWD/temp
        export PKGCACHE_OFFLINE=true

        create_riboseq_annotation.R \
            ${reference_twobit} \
            ${reference_gtf} \
            "\${PWD}/bsgenome_install" \
            "custom" \
            \${PWD} \
            ${reference_fasta}
        """
}

// Run RiboseQC on every sample
process riboseqc {

    tag "${sample_id}"
    label "Ribo_Seq_R"

    input:
        tuple val(sample_id), path(bam)
        path orfquant_annotation
        path package_install_loc
        val readlength_choice_method

    output:
        tuple val(sample_id), path("${sample_id}/${sample_id}_for_ORFquant"), emit: orfquant_psites
        path "${sample_id}/${sample_id}_results_RiboseQC_all", emit: riboseqc_all
        tuple path("${sample_id}/${sample_id}_P_sites_minus.bedgraph"), 
            path("${sample_id}/${sample_id}_P_sites_plus.bedgraph"),
            path("${sample_id}/${sample_id}_P_sites_uniq_minus.bedgraph"),
            path("${sample_id}/${sample_id}_P_sites_uniq_plus.bedgraph"), emit: bedgraphs
        path "${sample_id}/${sample_id}*"

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        run_riboseqc.R \
            ${bam} \
            ${sample_id}/${sample_id} \
            ${orfquant_annotation} \
            ${package_install_loc} \
            ${readlength_choice_method}
        """
}

// Sort RiboseQC output bedgraphs
process sort_bedgraphs{
    label "Ribo_Seq_tools"

    input:
        path bedgraph_file

    output:
        path "${bedgraph_file.simpleName}_sorted.bedgraph", emit: sorted_bedgraph

    script:
        """
        sort -k1,1 -k2,2n ${bedgraph_file} > ${bedgraph_file.simpleName}_sorted.bedgraph
        """
}

// Merge sorted bedgraphs into the correct groups
process merge_bedgraphs{

    label "Ribo_Seq_tools"

    input:
        path bedgraphs

    output:
        path "merged*.bedgraph", emit: merged_bedgraphs

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        declare -A groups
        groups[uniq_minus]="P_sites_uniq_minus"
        groups[uniq_plus]="P_sites_uniq_plus"
        groups[minus]="P_sites_minus"
        groups[plus]="P_sites_plus"

        # Iterate over every group
        for key in "\${!groups[@]}"; do
            pattern="\${groups[\$key]}"
            files=()
            # Check if file fits in group
            for file in ${bedgraphs}; do
                [[ "\$file" == *"\$pattern"* ]] && files+=( "\$file" )
            done
            echo \${files}
            # Merge the bedtools in each group
            bedtools unionbedg -i "\${files[@]}" > unsorted_merged_\${pattern}_sorted.bedgraph
            # Sum the p-site values into one big value
            awk '{sum=0; for (i=4; i<=NF; i++) sum+=\$i; print \$1"\t"\$2"\t"\$3"\t"sum}' unsorted_merged_\${pattern}_sorted.bedgraph > merged_\${pattern}.bedgraph
        done
        """
}

// Convert the bedgraph files to bigwig files
process convert_to_bigwig{

    label "Ribo_Seq_tools"

    input:
        path summed_bedgraph
        val genome_fai

    output:
        path "${summed_bedgraph.simpleName}.bw", emit: bigwig_p_site_tracks

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        # Create genome sizes file required for conversion to BigWig format
        cut -f1,2 ${genome_fai} > genome.sizes
        bedGraphToBigWig "${summed_bedgraph}" "genome.sizes" "${summed_bedgraph.simpleName}.bw" 
        """   
}

// Create RiboseQC html report
process create_riboseqc_report{

    label "Ribo_Seq_R"

    input:
        val riboseqc_all
        val html_template

    output:
        path "RiboseQC_report.html", emit: riboseqc_html_report

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    riboseqc_html.R \
    ${riboseqc_all.join(' ')} \
    ${html_template}
    """

}