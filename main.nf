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
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .set { ref_sdf_files }

Channel
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .set { ref_sdf_files2 }
    
Channel
    .fromPath("${params.pdb_sdf_files}/*.pdb")
    .set { pdb_files }

Channel
    .fromPath("${params.diffdock_location}/*", type: 'any')
    .set { diffd_tool }

Channel
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .map { [it.simpleName.split("__")[0], it.simpleName.split("__")[1].split("_")[0], it.simpleName] }
    .set { identifiers }

    
/*
* include the modules
*/

include { prepare_reference_files } from "./modules/prepare_reference_files"
include { prepare_ligand_sdf } from "./modules/prepare_ligand_sdf"
include { add_Hs_to_receptor } from "./modules/add_Hs_to_receptor"
include { p2rank } from "./modules/p2rank"
include { calculate_boxSize } from "./modules/calculate_boxSize"
include { docking_box } from "./modules/docking_box"
include { diffdock; diffdock_single; create_diffdock_csv } from "./modules/diffdock"
include { vina_prepare_receptor2; vina_prepare_ligand2; vina_box2; vina3; vina_pdbtqToSdf3 } from "./modules/vina"
include { rmsd } from "./modules/scoring"

/*
* main workflow
*/

workflow {
    
    /* 
    * preparation of reference and input for docking 
    */
    
    sdf_for_docking = prepare_ligand_sdf(ref_sdf_files.collect())
    
    pdb_files.filter{(it.simpleName =~ /_/)}.map{[it.baseName, it]}.set{ pdbs }
    pdb_Hs = add_Hs_to_receptor(pdbs)
   
    binding_pockets = p2rank(pdb_Hs)
    
    // define box parameters for vina-like tools
    sdf_for_docking.flatten().filter{!(it.simpleName =~ /_/)}.set { ligand_only }
    box_size = calculate_boxSize(ligand_only.map{file -> tuple(file.simpleName, file)})
    
    identifiers.combine(binding_pockets.pockets, by: 0)
                .combine(box_size, by: 1)
                .set{ input_dockingBox }
    boxes = docking_box(input_dockingBox)
   
   
    /* 
    * docking using Diffdock
    */
    
    diffd_csv = create_diffdock_csv(ref_sdf_files2.collect())
    diffdock_predictions = diffdock(diffd_csv, pdb_Hs.flatten().filter{it =~ /\//}.collect(), sdf_for_docking.collect(), diffd_tool.collect())
    
    //rmsd_out = rmsd(diffdock_predictions.predictions.collect())
    
    // diffdock single samples
    //ref_sdf_files.map{ [it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName] }
    //         .combine(pdb_Hs, by: 0)
    //         .combine(sdf_for_docking.flatten().map{file -> tuple(file, file.simpleName)}, by: 1)
    //         .set{ input_diffd_single }
    //diffdock_predictions = diffdock_single(input_diffd_single, diffd_tool.collect())
    //rmsd_out = rmsd(diffdock_predictions.predictions.collect())
    
    
    /* 
    * docking using Vina
    */
    
    preped_ligands = vina_prepare_ligand2(ligand_only.collect())
    preped_receptors = vina_prepare_receptor2(pdb_Hs)

    identifiers.combine(preped_receptors, by: 0)
            .combine(preped_ligands.flatten().map{file -> tuple(file, file.simpleName)}, by: 1)
            .map { [ it[2], it[0], it[1], it[3], it[4] ] }
            .combine(boxes.flatten().map{file -> tuple(file.simpleName.toString().split("_pocket")[0], file)}, by: 0)
            .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }
            .set {vina_input}

    vina_out = vina3(vina_input)
    vina_sdf = vina_pdbtqToSdf3(vina_out.vina_result)
    
}
