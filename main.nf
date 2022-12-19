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
    .value("${params.dataset}")
    .set { dataset }

Channel
    .value("${params.protein_ligand_csv}")
    .set { protein_ligand_csv }
   
//Channel
//    .value("$PWD/diffdock_predictions")
//    .set { diffdock_predictions }

Channel
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .map { file -> tuple(file.simpleName, file) }
    .set { sdf_files }

    
/*
* include the modules
*/

include { prepare_pdb_sdf } from "./modules/prepare_pdb_sdf"
include { diffdock; create_csv; diffdock_single } from "./modules/diffdock"
include { rmsd } from "./modules/scoring"

/*
* main workflow
*/

workflow {

    //pdb_sdf_files = prepare_pdb_sdf(dataset)
    diffdock_input_csv = create_csv(sdf_files)
    diffdock_predictions = diffdock_single(diffdock_input_csv)
    //diffdock_predictions = diffdock(pdb_sdf_files.protein_ligand_csv)
    //diffdock_predictions = diffdock(protein_ligand_csv)
    //rmsd_out = rmsd(diffdock_predictions)

}
