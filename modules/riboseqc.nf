// Run RiboseQC on every sample
process riboseqc {

    tag "${sample_id}"
    label "Ribo_Seq_R_scripts"
    publishDir "${outdir}/riboseqc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    val outdir
    val orfquant_annotation
    val package_install_loc

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_for_ORFquant"), emit: orfquant_psites
    path "${sample_id}/${sample_id}_results_RiboseQC_all", emit: riboseqc_all
    tuple path("${sample_id}/${sample_id}_P_sites_minus.bedgraph"), 
        path("${sample_id}/${sample_id}_P_sites_plus.bedgraph"),
        path("${sample_id}/${sample_id}_P_sites_uniq_minus.bedgraph"),
        path("${sample_id}/${sample_id}_P_sites_uniq_plus.bedgraph"), emit: bedgraphs
    path "${sample_id}/${sample_id}*"

    script:
    """
    run_riboseqc.R \
        ${bam} \
        ${sample_id}/${sample_id} \
        ${orfquant_annotation} \
        ${package_install_loc}
    """
}

// Sort RiboseQC output bedgraphs
process sort_bedgraphs{

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

    input:
    path bedgraphs

    output:
    path "merged*.bedgraph"

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
    // Publish to correct group ouput dir
    publishDir "${outdir}/igv_files/psite_tracks/minus", mode: 'copy', pattern: "*P_sites_minus*.bw"
    publishDir "${outdir}/igv_files/psite_tracks/plus", mode: 'copy', pattern: "*P_sites_plus*.bw"
    publishDir "${outdir}/igv_files/psite_tracks/uniq_minus", mode: 'copy', pattern: "*P_sites_uniq_minus*.bw"
    publishDir "${outdir}/igv_files/psite_tracks/uniq_plus", mode: 'copy', pattern: "*P_sites_uniq_plus*.bw"

    input:
    path summed_bedgraph
    val genome_fai
    val outdir

    output:
    path "${summed_bedgraph.simpleName}.bw"

    script:
    """
    # Create genome sizes file required for conversion to BigWig format
    cut -f1,2 ${genome_fai} > genome.sizes
    bedGraphToBigWig "${summed_bedgraph}" "genome.sizes" "${summed_bedgraph.simpleName}.bw" 
    """   
}
