# **Nextflow Ribo-seq Pipeline**

This pipeline is designed for the analysis and interpretation of Ribosome profiling data.

## **Requirements**
* Nextflow 23.04.1 or later
* Java 17 or later (up to 24)

### Containerised software
* Trimgalore         (0.6.6)
* Bowtie2            (2.5.4)
* STAR               (2.7.8)
* SAMtools           (1.12)
* ORFquant           (4.1.2)
* java (as a module) (1.8.0)
* bedgraphtobigwig   (ucsc 482)
* Bedtools           (2.31.0)
* TRISTAN            (1.0.0)
* MultiQC            (1.30)
* R, including the following packages:
    * tidyverse
    * tibble
    * dplyr
    * seqinr
    * stringdist
    * magrittr
    * GenomicRanges
    * rtracklayer
    * stringr
    * AnnotationDbi
    * Biostrings
    * GenomicFeatures
    * ggplot2
    * data.table
    * BSgenome
    * BSgenomeForge
    * txdbmaker
    * github repo = 'damhof/RiboseQC'
    * github repo = 'damhof/ORFquant

## **Overview**

1. **Quality Control, Trimming and Filtering**: Trim reads with Trimgalore and use a custom Bowtie2 index to remove RNA contaminants, keeping mRNA ribosome protected fragments.
2. **Alignment**: Use STAR to align the RPFs to the chosen reference genome and transcriptome.
3. **QC stats**: RiboseQC and multiQC for the output of relevant RPF statistics and the location of P-sites.
4. **ORF calling**: PRICE, ORFquant and RiboTIE algorithms for calling ORFs.
5. **ORF annotation**: Create annoated ORF table using ORF caller output.
6. **ORF harmonisation**: Merge the ORF annotation tables of the different ORF callers.
7. **ORF expression**: Obtain expression info of ORFs using the input samples, create a P-site and a PPM table.

## **Usage**

The pipeline can be run with a command as follows:

1. **Configure Parameters**:

```bash
nextflow run nf_riboseq_pipeline/main.nf \
    -c /path/to/local/params.config \
    -profile slurm
```

    Configure the parameters located in `params.config` ( which can be found in documentation) to suit your analysis.

    - Inputs: Adjust settings to toggle pipeline components on and off.

    - References and other files**: Define (file) paths to reference files, packages and other settings used in the pipeline.

    - Containers : Currently, we have loose containers for the entire pipeline that are described in `config/base.config`.


## Inputs

### Samplesheet

- The input samplesheet is a structured file, in CSV format, that lists all the samples to be processed. Each row represents a sample with detailed information required by the pipeline. Below is a breakdown of the expected columns and their constraints:

| **Column Name**  | **Description**                                                                 | **Required** | **Constraints**                                                                 |
|------------------|---------------------------------------------------------------------------------|--------------|---------------------------------------------------------------------------------|
| `subject_id`     | Unique identifier for the subject (e.g., patient).                              | Yes          | Must be a string without spaces.                                                |
| `sample_id`      | Unique identifier for the sample (e.g., sample barcode).                        | Yes          | Must be a string without spaces.                                                |
| `group_id`       | (Optional) Identifier for the sample group (e.g., cohort).                      | No           | Must be a string without spaces. If left empty, it must be omitted entirely.    |
| `sample_type`    | Type of sample: either `tumor` or `normal`.                                     | Yes          | Must be either `tumor` or `normal`.                                             |
| `sequence_type`  | Specifies the type of sequencing data: `rna` or `dna` or `ribo`.                | Yes          | Must be `ribo`, `rna` or `dna`.                                                  |
| `file_type`      | Format of the input files, e.g., `fastq`, `bam`, `cram`, `vcf`, `csv`, `tsv`.   | Yes          | Must be one of the supported formats: `fastq`, `bam`, `cram`, `vcf`, `csv`, etc.|
| `filename_1`     | Path to the first file (e.g., R1 FASTQ file for paired-end or single-end data). | Yes          | Must be a valid file path with no spaces. File extension must match `file_type`. |


**Example samplesheet**:
```
subject_id,sample_id,group_id,sample_type,sequence_type,file_type,filename_1
subject1,sampleA,cohort1,tumor,rna,fastq,/path/to/sampleA_R1.fastq.gz,/path/to/sampleA_R2.fastq.gz
subject2,sampleB,,normal,rna,bam,/path/to/sampleB.bam,
subject3,sampleC,cohort2,tumor,dna,vcf,/path/to/sampleC.vcf,
subject4,sampleA,cohort1,tumor,ribo,fastq,/path/to/sampleA.fastq.gz
```

Notes:
- All file paths (filename_1) must not contain spaces and should have extensions that match the declared file_type.
- Only samples that have ribo as sample_type, and fastq as file_type will be used in this pipeline.

### Parameter specification
### Input/output files

| Parameter      | Description                                                                                     |
| :--------------- | ------------------------------------------------------------------------------------------------- |
| --input        | **(Required)** Path to directory with raw FASTQ files (accepts a subdirectories one level deep) |
| --sample_sheet | **(Required)** CSV file with sample metadata (see format below)                                 |
| --outdir       | **(Required)** Path to output directory                                                         |

### Reference files


| Parameter              | Description                                           |
| :--------------------- | ----------------------------------------------------- |
| --reference_fasta      | **(Required)** FASTA file with reference genome       |
| --reference_fasta_fai  | **(Required)** FASTA file with reference genome fai   |
| --reference_twobit     | **(Required)** 2bit file with reference genome        |
| --reference_gtf        | **(Required)** GTF file with reference annotations    |
| --reference_protein_fa | **(Required)** FASTA file with reference proteins     |
| --contaminants_fasta   | **(Required)** FASTA file with known RNA contaminants |
| --package_install_loc  | **(Required)** Path to BSgenome directory             |
| --orfquant_annotation  | **(Required)** Path to ORFquant index .rannot file    |
| --bowtie2_index_path   | **(Optional)** Path to Bowtie2 index directory        |
| --price_index_path     | **(Optional)** Path to PRICE index file               |
| --star_index_path      | **(Optional)** Path to star_index directory           |

Note: the last three index paths can be automatically generated if the paths are empty

### Pipeline section Toggles


| Parameter    | Description                                  | Default |
| :------------------- | -------------------------------------------- | --------- |
| --run_qc             | Enable fastq trimming and RiboseQC           | true      |
| --run_orf_prediction | Enable ORF prediction and annotation         | true      |

### Pipeline module Toggles


| Parameter        | Description                                  | Default |
| :--------------- | -------------------------------------------- | ------- |
| --run_selection  | Run quality control and contaminants filter  | true    |
| --run_alignment  | Run STAR alignments                          | true    |
| --run_riboseqc   | Run RiboseQC on STAR output                  | true    |
| --run_orfquant   | Run ORFquant ORF calling                     | true    |
| --run_price      | Run PRICE ORF calling                        | true    |
| --run_ribotie    | Run RiboTIE (TRISTAN) ORF calling            | true    |
| --run_psite      | Run psite annotation                         | true    |
| --run_annotation | Run ORF annotation and harmonisation         | true    |
| --run_expression | Run ORF expression                           | true    |
| --run_multiqc    | Create MultiQC output report                 | true    |

Note: You can start the pipeline in a later step by setting these parameters. It will search the given outdir for the required files, make sure to give these files in the correct format and directory.

### Pipeline extra Toggles


| Parameter                  | Description                                                | Default |
| :------------------------- | ---------------------------------------------------------- | --------- |
| --help                     | Show pipeline options ( currently not working)             | false   |
| --validate_gtf_attributes  | Check first 10k gtf rows for correct format and attributes | false   |
| --keep_bam                 | Run RiboseQC on STAR output                                | true    |


### Optional Parameters

#### RiboTIE


| Parameter             | Description                                                                       | Default |
| :---------------------| --------------------------------------------------------------------------------- | ------- |
| --ribotie_min_samples | Minimal amount of samples RiboTIE needs to call the same ORF in for it to be kept | 2       |


## Outputs

The pipeline generates output files including trimmed reads, filtered reads, alignment results, RiboseQC output, ORF calls, ORF annotation tables, ORF expression tables, and a basic MultiQC report. The output directory is specified by the parameter `outdir` and has the following structure:

```
{outdir}
├── annotate_orfs
├── annotation
├── bedfiles
├── bowtie2
├── fastqc
├── final_orf_table
├── harmonise_orfs
├── igv_files
├── merged_ribotie
├── multiqc
├── orf_expression
├── orfquant
├── price
├── qc
├── riboseqc
├── ribotie
├── samplesheet
└── star
```

Additionally it produces Nextflow execution reports in `{project_folder}/log`

### Harmonised ORF table

The final_orf_table directory in the output directory holds the harmonised ORF table output, the following is a table explaining the columns:

|Column_name | Description|
| :--------------------- | -------------------------- |
|orf_id | The id of the ORF as given by the ORFcaller|
|summary_id | A combination of ORF location, ORF biotype, and its orf_id, to have a summary of relevant info in one column|
|gene_id | The gene_id based on the transcript id on which the ORFcaller says the ORF is on|
|gene_name | The gene_name based on the transcript id on which the ORFcaller says the ORF is on|
|gene_biotype | The gene_biotype based on the transcript id on which the ORFcaller says the ORF is on|
|protein_sequence | The protein sequence of the ORF|
|protein_length | The length of the protein sequence of the ORF|
|start_codon | The start codon of the ORF|
|stop_codon | The stop codon on the ORF|
|chr | The sequence on which the ORF is located|
|orf_start | The start coordinate of the ORF|
|orf_end | The stop coordinate of the ORF|
|strand | The strand of the ORF|
|starts | All the starts of the ORF CDS seperated by "_"|
|ends | All the ends of the ORF CDS seperated by "_"|
|has_p0_cds_overlap | Boolean, set to true if the ORF has inframe overlap with a reference CDS. Note that the ORF doesn't have to fit the exon boundaries of the transcript this CDS comes from to be set to TRUE.|
|cds_isoform_transcripts | All the transcripts of which the ORF could be an isoform off. This means that the ORF does fully fit within the exon boundaries of this transcript, but the exon length from the start of the ORF to the stop of the ORF is longer than the ORF length|
|orfcaller | The used ORFcaller (currently ORFquant, PRICE, or RiboTIE)|
|tx_id | All transcripts on which this ORF can be located, meaning it fits within the exon boundaries of the given transcripts, transcript ids are seperated by "__"|
|transcript_biotype_all | The transcript_biotypes of all transcripts in the tx_id column, transcript biotypes seperated by "__"|
|orf_biotypes_all | The ORF biotype (category) of the ORF classified on the transcripts in tx_id, shown in the same order as the transcripts in tx_id, seperated by "__"|
|orf_biotype_single | The single ORF biotype (category) chosen from orf_biotypes_all. This single biotype is chosen using orf_biotypes_all and a factor in which it will pick the first occurence|

**Further orf_biotype_single explanation**
The orf_biotype_single is decided by comparing the orf_biotypes_all to the factor shown below, in which it will pick the first occurence of the factor that is found in orf_biotypes all. This means e.g. that even if there is only one transcript which is classified as an ORF-annotated in orf_biotypes_all the orf_biotype_single will become ORF-annotated even if there are other transcripts on which the ORF would have a different ORF biotype. 

**The used factor:**
```
"ORF-annotated", "NC-variant", 
"ovCDS_lncRNA-ORF", "lncRNA-ORF", 
"pseudogene-ORF", "Incomplete_overlap", 
"intORF", "uoORF", "doORF", "uORF", "dORF",
"Processed_transcript_ORF", "novel-ORF",
"undefined"
```

### The ORF biotypes and their explanation:

|ORF biotype | Explanation|
| :--------- | ---------- |
|ORF-annotated | An ORF which matches the structure of a canonical protein|
|NC-variant | An ORF which overlaps a canonical protein but has a N or C terminal truncation or extension. Can have a combination of these options|
|ovCDS_lncRNA-ORF | An ORF found on a transcript which has transcript type lncRNA, but it does have inframe overlap with a CDS|
|lncRNA-ORF | An ORF found on a transcript which has transcript biotype lncRNA|
|pseudogene-ORF | An ORF found on a transcript which has transcript biotype pseudogene|
|intORF | An ORF found within the exon boundaries of a transcript but not in the same frame|
|uoORF | An upstream overlapping ORF which is not inframe with reference CDS|
|doORF | A downstream overlapping ORF which is not inframe with reference CDS|
|uORF | An upstream ORF which is not inframe with reference CDS|
|dORF | A downstream ORF which is not inframe with reference CDS|
|Processed_transcript_ORF | An ORF found within the exon boundaries of a transcript for which no CDS is found|
|novel-ORF | An ORF of which the transcript was created by StringTie|
|undefined | An ORF without transcript assignment|



## Support and Contributions

- For issues or questions, [create an issue](https://github.com/VanHeeschTools/nf_ribo_pipeline/issues).
