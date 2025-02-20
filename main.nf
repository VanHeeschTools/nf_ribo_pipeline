include { paramsHelp } from 'plugin/nf-schema'
include { printHeader } from "./modules/helperFunctions.nf"
include { RIBOSEQ } from "./workflows/RIBOSEQ.nf"

workflow {

    printHeader()

    if (params.help) {

        log.info paramsHelp("nextflow run VanHeeschTools/nf_ribo_pipeline -c params.config")

    } else {

        RIBOSEQ()
    }

}
