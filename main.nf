#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Define the pipeline parameters
 */

log.info """
DTBW  ~  version ${workflow.manifest.version}
=============================================
input directory        : ${params.pdb_sdf_files}
input naming           : ${params.naming}
receptor_Hs            : ${params.receptor_Hs}
diffdock_mode          : ${params.diffdock_mode}
"""


/*
* define the input channels
*/

Channel
    .fromPath("${params.pdb_sdf_files}/*.sdf")
    .set { ref_sdf_files }

Channel
    .fromPath("${params.pdb_sdf_files}/*.mol2")
    .set { mol_files }

Channel
    .fromPath("${params.pdb_sdf_files}/*.pdb")
    .set { pdb_files }

Channel
    .fromPath("${params.diffdock_tool}/*", type: 'any')
    .set { diffd_tool }

if (params.receptor_Hs == "yes") {
    Channel
        .fromPath("${params.pdb_sdf_files}/*.pdb")
        .map { [it.simpleName.split("_")[0], it] }
        .set { pdb_Hs }
}


/*
* define identifiers to combine the files: receptor name , ligand name, complex name
*/

if (params.naming == "default") {
    ref_sdf_files.map { [it.simpleName.split("__")[0], it.simpleName.split("__")[1].split("_")[0], it.simpleName] }
                 .set { identifiers }
}
else {
    ref_sdf_files.map { [it.simpleName.split("_")[0], it.simpleName, it.simpleName] }
                 .set { identifiers }
}


/*
* include the modules
*/

include { prepare_ligand_sdf } from "./modules/prepare_ligand_sdf"
include { add_Hs_to_receptor } from "./modules/add_Hs_to_receptor"
include { p2rank } from "./modules/p2rank"
include { calculate_boxSize } from "./modules/calculate_boxSize"
include { docking_box } from "./modules/docking_box"
include { create_diffdock_csv; diffdock; diffdock_single } from "./modules/diffdock"
include { vina_prepare_receptor; vina_prepare_ligand; vina; pdbtqToSdf as vina_pdbtqToSdf; pdbtqToSdf as smina_pdbtqToSdf; pdbtqToSdf as gnina_pdbtqToSdf } from "./modules/vina"
include { gnina } from "./modules/gnina"
include { smina } from "./modules/smina"
include { tankbind } from "./modules/tankbind"
include { ost_scoring as tb_ost; ost_scoring as dd_ost; ost_scoring as vina_ost; ost_scoring as smina_ost; ost_scoring as gnina_ost } from "./modules/scoring"


/*
* main workflow
*/

workflow {

    /*
    * preprocessing
    */

    sdf_for_docking = prepare_ligand_sdf(ref_sdf_files.collect(), mol_files.collect().ifEmpty([]))

    if (params.naming == "default") {
        sdf_for_docking.sdf_files.flatten().filter{!(it.simpleName =~ /_/)}.set { ligand_only }
        ligand_only.map{ [it.simpleName, it] }.set { ligand_tuple }
    }
    else {
        sdf_for_docking.sdf_files.flatten().set{ ligand_only }
        ligand_only.map{ [it.simpleName.toString().split("_preped")[0], it] }.set{ ligand_tuple }
    }

    if (params.receptor_Hs == "no") {
        pdb_files.filter{(it.simpleName =~ /_/)}.map{[it.baseName, it]}.set{ pdbs }
        pdb_Hs = add_Hs_to_receptor(pdbs)
    }

    /*
    * binding pocket prediction
    */

    binding_pockets = p2rank(pdb_Hs)

    // define box parameters for vina-like tools

    box_size = calculate_boxSize(ligand_tuple)

    identifiers.combine(pdb_Hs, by: 0)
               .combine(binding_pockets.pockets, by: 0)
               .combine(box_size.size, by: 1)
               .set{ input_dockingBox }
    boxes = docking_box(input_dockingBox)

    /*
    * docking using Diffdock
    */

    if (params.diffdock_mode == "batch") {
        diffd_csv = create_diffdock_csv(ref_sdf_files.collect())
        diffdock_predictions = diffdock(diffd_csv, pdb_Hs.flatten().filter{it =~ /\//}.collect(), sdf_for_docking.sdf_files.collect(), diffd_tool.collect())
    }
    else if (params.diffdock_mode == "single") {
        if (params.naming == "default") {
            ref_sdf_files.map{ [it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName] }
                     .combine(pdb_Hs, by: 0)
                     .combine(sdf_for_docking.sdf_files.flatten().map{file -> tuple(file, file.simpleName)}, by: 1)
                     .set{ input_diffd_single }
        }
        else {
            identifiers.combine(pdb_Hs, by: 0)
                   .combine(ligand_tuple.map{ [ it[1], it[0] ] }, by: 1)
                   .set{ input_diffd_single }
        }

        diffdock_predictions = diffdock_single(input_diffd_single, diffd_tool.collect())
    }

    /*
    * input preparation vina-like tools
    */

    preped_ligands = vina_prepare_ligand(ligand_tuple)
    preped_receptors = vina_prepare_receptor(pdb_Hs)

    identifiers.combine(preped_receptors, by: 0)
               .combine(preped_ligands.map{ [it[1], it[0]] }, by: 1)
               .map { [ it[2], it[0], it[1], it[3], it[4] ] }
               .combine(boxes.flatten().map{file -> tuple(file.simpleName.toString().split("_pocket")[0], file)}, by: 0)
               .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }
               .set {vina_input}

    /*
    * docking using Vina
    */

    vina_out = vina(vina_input)
    vina_sdf = vina_pdbtqToSdf(vina_out.vina_result, Channel.value( 'vina' ))

    /*
    * docking using smina
    */

    smina_out = smina(vina_input)
    smina_sdf = smina_pdbtqToSdf(smina_out.smina_result, Channel.value( 'smina' ))

    /*
    * docking using gnina
    */

    gnina_out = gnina(vina_input)
    gnina_sdf = gnina_pdbtqToSdf(gnina_out.gnina_result, Channel.value( 'gnina' ))

    /*
    * docking using tankbind
    */

    identifiers.combine(pdb_Hs, by: 0)
               .combine(binding_pockets.pockets, by: 0)
               .combine(ligand_only.map{ [it, it.simpleName.toString().split("_preped")[0]] }, by: 1)
               .set {tankbind_input}

    tankbind_out = tankbind(tankbind_input)


    /*
    * ost scoring
    */

    if (params.naming == "default") {
        pdb_files.map{ [it.baseName, it] }.set{ ref_pdbs }
    }
    else {
        pdb_files.map{ [it.simpleName.split("_")[0], it] }.set{ ref_pdbs }
    }

    identifiers.combine(ref_pdbs, by: 0).map { [ it[2], it[0], it[1], it[3] ] }         // complex, receptor, ligand, rec_pdb file
               .combine(ref_sdf_files.map { [it.simpleName, it] }, by: 0)    // complex, receptor, ligand, rec_pdb file, ref_lig file
               .set { reference_files }

    // tankbind

    reference_files.combine(tankbind_out.sdfs, by: 0)
                   .set { tb_scoring_input }
    tb_scores = tb_ost(tb_scoring_input, Channel.value( 'tankbind' ))

    // diffdock
    diffdock_predictions.predictions.flatten().map{[it.toString().split('/')[-2], it]}.groupTuple().view()
    reference_files.combine(diffdock_predictions.predictions.flatten().map{[it.toString().split('/')[-2], it]}.groupTuple(), by: 0)
                    .set { dd_scoring_input }
    dd_scores = dd_ost(dd_scoring_input, Channel.value( 'diffdock' ))

    // vina
    reference_files.combine(vina_sdf.map{[ it[0], it[3]]}, by: 0)
                   .set { vina_scoring_input }
    vina_scores = vina_ost(vina_scoring_input, Channel.value( 'vina' ))


    // smina
    reference_files.combine(smina_sdf.map{[ it[0], it[3]]}, by: 0)
                   .set { smina_scoring_input }
    smina_scores = smina_ost(smina_scoring_input, Channel.value( 'smina' ))


    // gnina
    reference_files.combine(gnina_sdf.map{[ it[0], it[3]]}, by: 0)
                   .set { gnina_scoring_input }
    gnina_scores = gnina_ost(gnina_scoring_input, Channel.value( 'gnina' ))
}
