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
    .map { [it.baseName.toString().split("__")[0], it.baseName.toString().split("__")[1].split("_")[0], it.simpleName] }
    .set { identifiers }

Channel
    .fromPath("${params.pdb_sdf_files}/*.pdb")
    .map { file -> tuple(file.simpleName, file) }
    .filter{!(it[0] =~ /_/)}
    .flatten().filter{it =~ /\//}
    .set { molecules }
    

// above tested ok
//
//Channel
//    .value("${params.protein_ligand_csv}")
//    .set { protein_ligand_csv }
   
//Channel
//    .value("$PWD/diffdock_predictions")
//    .set { diffdock_predictions }
    
//Channel
//    .fromPath("${params.pdb_sdf_files}/*.sdf")
//    .map { file -> tuple(file.simpleName, file) }
//    .set { sdf_files }

Channel
    .fromPath("${params.pdb_sdf_files}/")
    .set { pdb_sdf_files }
    
Channel
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .map { [it.simpleName, it.baseName.toString().split("__")[0], it] }
    .set { ligands }
    
Channel
    .value("${params.pdb_sdf_files}/")
    //.fromPath("/scicore/home/schwede/leeman0000/projects/LIGATE/test_vina/data_tutorial/*.pdb")
    //.map { file -> tuple(file.simpleName, file) }
    .set { receptors }

    
/*
* include the modules
*/

include { prepare_reference_files } from "./modules/prepare_reference_files"
include { prepare_ligand_sdf } from "./modules/prepare_ligand_sdf"
include { diffdock; diffdock_single; create_diffdock_csv } from "./modules/diffdock"
include { rmsd } from "./modules/scoring"
include { vina_prepare_receptor2; vina_prepare_ligand2; vina_box2; vina2; vina_pdbtqToSdf } from "./modules/vina"

/*
* main workflow
*/

workflow {
    
    /* 
    * preparation of reference and input for docking 
    */
    
    //ref_pdb_sdf_files = prepare_reference_files(dataset)
    sdf_for_docking = prepare_ligand_sdf(ref_sdf_files.collect())
   
   
    /* 
    * docking using Diffdock
    */
    
    diffd_csv = create_diffdock_csv(ref_sdf_files2.collect())
    diffdock_predictions = diffdock(diffd_csv, pdb_files.collect(), sdf_for_docking.collect(), diffd_tool.collect())
    
    //rmsd_out = rmsd(diffdock_predictions)
    
    // diffdock single samples
    //diffdock_predictions = diffdock_single(sdf_files, diffd_tool.collect())
    //rmsd_out = rmsd(diffdock_predictions.predictions.collect())
    
    
    /* 
    * docking using Vina
    */
    
    sdf_for_docking.flatten().filter{!(it.simpleName =~ /_/)}.set { ligand_only }
    preped_ligands = vina_prepare_ligand2(ligand_only.collect())
    preped_receptors = vina_prepare_receptor2(pdb_files.filter{(it.simpleName =~ /_/)}.collect())
    vina_box_out = vina_box2(pdb_files.filter{(it.simpleName =~ /_/)}.map{file -> tuple(file.simpleName, file)}, pdb_files.filter{!(it.simpleName =~ /_/)}.collect())

    identifiers.combine(preped_receptors.flatten().map{file -> tuple(file.simpleName, file)}, by: 0)
            .combine(vina_box_out.flatten().map{file -> tuple(file.simpleName.split("_b")[0], file)}, by: 0)
            .combine(preped_ligands.flatten().map{file -> tuple(file, file.simpleName)}, by: 1)
            .set {vina_input}
    
    vina_out = vina2(vina_input)
    vina_sdf = vina_pdbtqToSdf(vina_out.vina_result)
    
}
