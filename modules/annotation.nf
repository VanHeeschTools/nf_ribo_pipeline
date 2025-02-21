process annotate_orfs {

    label "annotate_orfs"
    publishDir "${outdir}/annotate_orfs", mode: 'copy'

    input:

    output:

    script:
    """
    orf_annotate.R \
    ""\
    ""


    #args <- commandArgs(trailingOnly = TRUE)
    #orfs_loc <- args[1]
    #txdb_loc <- args[2]
    #orfcaller <- args[3]
    #annotation_provider <- args[4]
    #gencode_uniprot_file <- args[5]
    #uniprot_protein_fasta_loc <- args[6]
    #cpus <- args[7]

    """
}

process harmonise_orfs{

    // Should only be run when both ORFquant and PRICE are run

    label "harmonise_ofs"
    publishDir "${outdir}/harmonise_orfs", mode: 'copy'

    input:

    output:

    script:
    """
    orf_harmonisation.R \
    ""\
    ""
    # Needs a gtf, txdb
    # Needs the annotation r script output for price and orfquant

    """

}