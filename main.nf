#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Define the pipeline parameters
 */

log.info """
PickyBinder  ~  version ${workflow.manifest.version}
=============================================
input directory        : ${params.pdb_sdf_files}
input naming           : ${params.naming}
AlphaFold models       : ${params.alphafold}
reference files        : ${params.ref_files}
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
    .fromPath("${params.ref_files}/*.{pdb,cif}")
    .set { scoring_ref }

Channel
    .fromPath("${params.diffdock_tool}/*", type: 'any')
    .set { diffd_tool }

if (params.receptor_Hs == "yes" && params.naming == "default") {
    Channel
        .fromPath("${params.pdb_sdf_files}/*.pdb")
        .map { [it.simpleName, it] }
        .set { pdb_Hs }
}
else if (params.receptor_Hs == "yes" && params.naming == "other") {
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

include { ligand_preprocessing } from "./modules/ligand_preprocessing"
include { add_Hs_to_receptor } from "./modules/add_Hs_to_receptor"
include { p2rank } from "./modules/p2rank"
include { calculate_boxSize } from "./modules/calculate_boxSize"
include { docking_box } from "./modules/docking_box"
include { create_diffdock_csv; diffdock; diffdock_single } from "./modules/diffdock"
include { vina_prepare_receptor; vina_prepare_ligand; vina; pdbtqToSdf as vina_pdbtqToSdf; pdbtqToSdf as smina_pdbtqToSdf; pdbtqToSdf as gnina_pdbtqToSdf } from "./modules/vina"
include { gnina_sdf } from "./modules/gnina"
include { smina_sdf } from "./modules/smina"
include { tankbind } from "./modules/tankbind"
include { ost_scoring as tb_ost; ost_scoring as dd_ost; ost_scoring as vina_ost; ost_scoring as smina_ost; ost_scoring as gnina_ost; score_summary } from "./modules/scoring"
include { ost_scoring_modelReceptors; combine_modelReceptors_scores } from "./modules/scoring"


/*
* main workflow
*/

workflow {

    /*
    * preprocessing
    */

    sdf_for_docking = ligand_preprocessing(ref_sdf_files.collect(), mol_files.collect().ifEmpty([]))

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
    if (params.tools =~ /vina/ || params.tools =~ /smina/ || params.tools =~ /gnina/) {
        box_size = calculate_boxSize(ligand_tuple)

        identifiers.combine(pdb_Hs, by: 0)
                   .combine(binding_pockets.pockets, by: 0)
                   .combine(box_size.size, by: 1)
                   .set{ input_dockingBox }
        boxes = docking_box(input_dockingBox)
    }


    /*
    * docking using Diffdock
    */

    if (params.tools =~ /diffdock/) {

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
    }
    else {
        dd_scores_for_summary = Channel.empty()
    }


    /*
    * docking using Vina
    */

    if (params.tools =~ /vina/) {
        preped_ligands = vina_prepare_ligand(ligand_tuple)
        preped_receptors = vina_prepare_receptor(pdb_Hs)

        identifiers.combine(preped_receptors, by: 0)
               .combine(preped_ligands.map{ [it[1], it[0]] }, by: 1)
               .map { [ it[2], it[0], it[1], it[3], it[4] ] }
               .combine(boxes.flatten().map{file -> tuple(file.simpleName.toString().split("_pocket")[0], file)}, by: 0)
               .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }
               .set {vina_input}

        vina_out = vina(vina_input)
        vina_sdf = vina_pdbtqToSdf(vina_out.vina_result, Channel.value( 'vina' ))
    }
    else {
        vina_scores_for_summary = Channel.empty()
    }


    /*
    * docking using smina
    */

    if (params.tools =~ /smina/ || params.tools =~ /gnina/) {
        identifiers.combine(pdb_Hs, by: 0)
                   .combine(ligand_tuple.map{ [it[1], it[0]] }, by: 1)
                   .map { [ it[2], it[0], it[1], it[3], it[4] ] }
                   .combine(boxes.flatten().map{file -> tuple(file.simpleName.toString().split("_pocket")[0], file)}, by: 0)
                   .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }
                   .set {smi_gni_input}
    }

    if (params.tools =~ /smina/) {
        smina_out = smina_sdf(smi_gni_input)
    }
    else {
        smina_scores_for_summary = Channel.empty()
    }


    /*
    * docking using gnina
    */

    if (params.tools =~ /gnina/) {
        gnina_out = gnina_sdf(smi_gni_input)
    }
    else {
        gnina_scores_for_summary = Channel.empty()
    }


    /*
    * docking using tankbind
    */

    if (params.tools =~ /tankbind/) {
        identifiers.combine(pdb_Hs, by: 0)
               .combine(binding_pockets.pockets, by: 0)
               .combine(ligand_only.map{ [it, it.simpleName.toString().split("_preped")[0]] }, by: 1)
               .set {tankbind_input}

        tankbind_out = tankbind(tankbind_input)
    }
    else {
        tb_scores_for_summary = Channel.empty()
    }


    /*
    * ost scoring
    */

    if (params.naming == "default") {
            scoring_ref.map{ [it.simpleName.split("_")[0], it] }
                         .combine(pdb_files.map{ [it.simpleName.split("_")[0], it.baseName, it] }, by: 0)
                         .map { [ it[2], it[1], it[3] ] }
                         .set{ ref_for_scoring }
    }
    else {
            scoring_ref.map{ [it.simpleName.split("_")[0], it] }
                         .combine(pdb_files.map{ [it.simpleName.split("_")[0], it] }, by: 0)
                         .set{ ref_for_scoring }
    }

    if (params.alphafold == "yes") {
            modelled_receptors_scores = ost_scoring_modelReceptors(ref_for_scoring)
            modelled_receptors_scores_summary = combine_modelReceptors_scores(modelled_receptors_scores.toList().flatten().filter{ it =~ /json/ }.collect())
    }

    identifiers.combine(ref_for_scoring, by: 0).map { [ it[2], it[0], it[1], it[3], it[4] ] }       // complex, receptor, ligand, ref pdb file, model pdb file
               .combine(ref_sdf_files.map { [it.simpleName, it] }, by: 0)    // complex, receptor, ligand, ref pdb file, model pdb file, ref_lig file
               .set { reference_files }

    // tankbind

    if (params.tools =~ /tankbind/) {
        reference_files.combine(tankbind_out.sdfs, by: 0)
                       .set { tb_scoring_input }
        tb_scores = tb_ost(tb_scoring_input, Channel.value( 'tankbind' ))
        tb_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ tb_scores_for_summary }
    }


    // diffdock
    if (params.tools =~ /diffdock/) {
        reference_files.combine(diffdock_predictions.predictions.flatten().map{[it.toString().split('/')[-2], it]}.groupTuple(), by: 0)
                       .set { dd_scoring_input }

        dd_scores = dd_ost(dd_scoring_input, Channel.value( 'diffdock' ))
        dd_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ dd_scores_for_summary }
    }


    // vina
    if (params.tools =~ /vina/) {
        reference_files.combine(vina_sdf.map{[ it[0], it[3]]}, by: 0)
                       .set { vina_scoring_input }
        vina_scores = vina_ost(vina_scoring_input, Channel.value( 'vina' ))
        vina_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ vina_scores_for_summary }
    }


    // smina
    if (params.tools =~ /smina/) {
        reference_files.combine(smina_out.smina_sdf.map{[ it[0], it[3]]}, by: 0)
                       .set { smina_scoring_input }
        smina_scores = smina_ost(smina_scoring_input, Channel.value( 'smina' ))
        smina_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ smina_scores_for_summary }
    }


    // gnina
    if (params.tools =~ /gnina/) {
        reference_files.combine(gnina_out.gnina_sdf.map{[ it[0], it[3]]}, by: 0)
                       .set { gnina_scoring_input }
        gnina_scores = gnina_ost(gnina_scoring_input, Channel.value( 'gnina' ))

        gnina_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ gnina_scores_for_summary }
    }

    // create score summary file
    overall_scores = score_summary(tb_scores_for_summary.ifEmpty([]).combine(
                                   dd_scores_for_summary.ifEmpty([]).combine(
                                   vina_scores_for_summary.ifEmpty([]).combine(
                                   smina_scores_for_summary.ifEmpty([]).combine(
                                   gnina_scores_for_summary.ifEmpty([]))))))
}
