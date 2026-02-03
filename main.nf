include { paramsHelp } from 'plugin/nf-schema'
include { printHeader } from "./modules/helperFunctions.nf"
include { RIBOSEQ } from "./workflows/RIBOSEQ.nf"

workflow {

    printHeader()
    
    RIBOSEQ()

}
