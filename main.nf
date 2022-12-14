#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Define the pipeline parameters
 */


// this prints the input parameters
log.info """
DTBW  ~  version ${workflow.manifest.version}
=============================================
"""

/*
* check the input type and create relevant channels
*/

Channel
    .value("${params.dataset})
    .set { dataset }

/*
* include the modules
*/

include { prepare_pdb_sdf } from "./modules/prepare_pdb_sdf"


/*
* main workflow
*/

workflow {

pdb_sdf_files = prepare_pdb_sdf(dataset)

}
