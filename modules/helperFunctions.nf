def printHeader () {
    def logMessage =  """
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        """
        log.info logMessage.stripIndent()
}

//TODO: make this work for the riboseq pipeline
def check_files(name, path, type) {
    if(!path) {
        error " File${name} of type ${type} not found at ${path}."
    } else {
        def file_to_check = null
        if (type == "dir") {
            file_to_check = file(path, type: "dir")
        } else {
            file_to_check = file(path, type: "file")
        }
        if (file_to_check != null) {
            // Check if it's a list of files
            if (file_to_check instanceof List) {
                // Check each file in the list
                file_to_check.flatten().each { file ->
                    if (!file.exists()) {
                        error "--${name}: ${type} doesn't exist, check path ${file}"
                    }
                }
            } else {
                // Check the single file
                if (!file_to_check.exists()) {
                    error "--${name}: ${type} doesn't exist, check path ${path}"
                }
            }
        } else {
            error "--${name}: No files found at ${path}"
        }
    }
}

def checkInputFiles() {
    // Locate bams
    def default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    def bam_avail = true
    if (!params.align) {
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}".stripIndent()
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}".stripIndent()
        } else {
            bam_avail = null
        }
    }

    // Check if container folder exists
    check_files("container_folder", params.container_folder, "file")

    // Check reference files
    check_files("reference_gtf", params.reference_gtf, "file")
    check_files("reference_fasta", params.reference_fasta, "file")
    check_files("reference_fasta_fai", params.reference_fasta_fai, "file")
    check_files("reference_protein_fa", params.reference_protein_fa, "file")
    check_files("twobit", params.twobit, "file")

    // Check if bowtie2 contaminants exist
    check_files("contaminants_fasta", params.contaminants_fasta, "file")

    // Check if ORFquant input files and directories exist
    check_files("package_install_loc", params.package_install_loc, "dir")
    check_files("orfquant_annotation", params.orfquant_annotation, "file")
    check_files("orfquant_annot_package", params.orfquant_annot_package, "dir")

    log.info "\n==========================\n"
}


// Check if all bowtie2 index paths exist and return its path if TRUE otherwise return null
def validate_bowtie2_index (String bowtie2_index){
    // Define the Bowtie2 index file extensions
    def bowtie2_extensions = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

    return bowtie2_extensions.every { ext ->
        file("${bowtie2_index}${ext}").exists()
    }
}

// Check if all STAR index paths exist and return its path if TRUE otherwise return null
def validate_star_index(String star_index_path){
    def star_index_files = [
        'chrLength.txt',
        'chrStart.txt',
        'geneInfo.tab',
        'SA',
        'sjdbList.fromGTF.out.tab',
        'chrNameLength.txt',
        'exonGeTrInfo.tab',
        'Genome',
        'SAindex',
        'sjdbList.out.tab',
        'chrName.txt',
        'exonInfo.tab',
        'genomeParameters.txt',
        'sjdbInfo.txt',
        'transcriptInfo.tab',
    ]

    // Check existence of all STAR index files
    return star_index_files.every { filename ->
            file("${star_index_path}/${filename}").exists()
    }
}

// Check if the PRICE index exists 
def validate_price_index(String price_index_path){
    return file("${price_index_path}/PRICE_index.oml").exists()
}

def validateGTF(String gtfPath) {
    def requiredAttrs = [
        gene      : ["gene_id", "gene_biotype", "gene_name"],
        transcript: ["gene_id", "gene_biotype", "gene_name", "transcript_id"],
        exon      : ["gene_id", "gene_biotype", "gene_name", "transcript_id", "exon_number"]
    ]

    def exonByTranscript = [:].withDefault { [] }
    def missingAttrs = []
    def seenTranscripts = new LinkedHashSet()

    println "Validating input GTF: ${gtfPath}"

    // Parse gtf and obtain attributes
    new File(gtfPath).eachLine { line ->
        if (line.startsWith("#") || !line.trim()) return

        def fields = line.split("\t")
        if (fields.size() < 9) return

        def type   = fields[2]
        def start  = fields[3].toInteger()
        def strand = fields[6]
        def attrText = fields[8]

        // Parse attributes into map
        def attrs = [:]
        attrText.split(";").each { part ->
            part = part.trim()
            if (!part) return
            def kv = part.split(/\s+/, 2)
            if (kv.size() == 2) {
                def key = kv[0]
                def value = kv[1].replaceAll(/^"|"$/, "")
                attrs[key] = value
            }
        }

        // Check if required attributes are present
        if (requiredAttrs.containsKey(type)) {
            def missing = requiredAttrs[type].findAll { !attrs.containsKey(it) || !attrs[it] }
            if (missing) {
                missingAttrs << "Missing ${missing.join(", ")} in ${type} line: ${line.take(80)}..."
            }
        }

        // Collect exons for first 10k transcripts
        if (type == "exon" && attrs.transcript_id) {
            def tid = attrs.transcript_id
            seenTranscripts << tid
            if (seenTranscripts.size() <= 10_000) {
                exonByTranscript[tid] << [start: start, strand: strand]
            }
        }
    }

    // Throw error if attributes are missing
    if (missingAttrs) {
        println "Attribute check failed:"
        missingAttrs.take(10).each { println "- $it" }
        if (missingAttrs.size() > 10)
            println "... and ${missingAttrs.size() - 10} more."
        System.exit(1)
    } else {
        println "All required attributes are present."
    }

    // Check if exons are sorted correctly
    exonByTranscript.each { tid, exons ->
        def strand = exons[0].strand
        def starts = exons*.start

        // Check that the starts follow correct strand-specific order
        def isOrdered = (strand == "+") ? 
            starts.inject([true, starts[0]]) { acc, val -> [acc[0] && val >= acc[1], val] }[0] :
            starts.inject([true, starts[0]]) { acc, val -> [acc[0] && val <= acc[1], val] }[0]
        // Throw error if not sorted correctly
        if (!isOrdered) {
            println "Exons not sorted correctly for transcript ${tid} on strand '${strand}'"
            println "Start positions: ${starts}"
            System.exit(1)
        }
    }
    println "Exons ordered correctly in first 10,000 transcripts."
    println "GTF validation passed for ${exonByTranscript.size()} transcripts."
}

// Copy an input samplesheet to outdir/samplesheet/
def copy_samplesheet(String input, String outdir) {
    // Check if input samplesheet and output directory paths are given
    if (!input || !outdir)
        return false

    // Check if input samplesheet file exists
    def inputFile = file(input)
    if (!inputFile.exists())
        return false

    // Creates the directory if needed
    def destDir = file("${outdir}/samplesheet")
    destDir.mkdirs()

    // Copy the samplesheet to output directory
    inputFile.copyTo(destDir.resolve(inputFile.name))

    return true
}

// Search output dir to collect ouput data from previous runs
def collect_output_previous_run(
    pattern,                    
    mode = 'sample_id',         
    emptyInsteadOfNull = false, 
    step = null                 
    ) {

    // Pattern: String, given path pattern to search for output files
    // Mode: String, either sample_id or path, should the output be a tuple of sample_id and ouput path or just the output path alone
    // emptyInsteadOfNull: Bool, should the output be an empty channel or null
    // step: String, log information about the step to write.

    // Check if given pattern gives any output files
    if ( file(pattern).isEmpty() ) {
        if (emptyInsteadOfNull){
            log.info "WARNING: No output found for ${step}, the output of this step will not be included in the rest of the pipeline."
        } else {
            log.info "ERROR: No matching output files found for step: ${step}, please make sure to provide them correctly or re-run the step."
        }
        // Return empty channel or null
        return emptyInsteadOfNull ? Channel.empty() : null
    }

    // If sample_id is required just create tuple channel of sample_ids and found files
    if (mode == 'sample_id') {
        return Channel.fromFilePairs(pattern, size: 1, checkIfExists: true)
    // If sample_id is not required create channel of found files
    } else if (mode == "path") {
        return Channel.fromPath(pattern)
    } else{
        return null 
    }
}