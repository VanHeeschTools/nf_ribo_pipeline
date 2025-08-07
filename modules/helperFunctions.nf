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
        error "When running merge without assembly you must provide `${name}`."
    } else {
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
    //Check inputs

    // Locate bams
    default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    bam_avail = true
    if (!params.align) {
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}".stripIndent()
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}".stripIndent()
        } else {
            bam_avail = null
        }
    }

    // Check reference_gtf
    check_files("reference_gtf", params.reference_gtf, "file")

    // Check contaminants fasta

    // Check bowtie2 index
    if (params.qc){
        check_files("kallisto_index", params.kallisto_index, "file")
    }

    // Check references for build_annotaiton
    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

    log.info "\n==========================\n"
}

process write_collected_paths {
    input:
    val collected_paths

    output:
    path "file_paths.txt", emit: collected_file_channel

    script:
    """
    printf "%s\n" "${collected_paths.join('\n')}" > file_paths.txt
    """
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
