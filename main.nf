#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Define the pipeline parameters
 */

log.info """
PickyBinder  ~  version ${workflow.manifest.version}
=============================================
input data             : ${params.data}
input naming           : ${params.naming}
reference files        : ${params.ref_files}
AlphaFold models       : ${params.alphafold}
receptor_Hs            : ${params.receptor_Hs}
diffdock_mode          : ${params.diffdock_mode}
"""


/*
* define the input channels
*/

// check if there is a path to the input data
if (params.data == "") {
    println()
    println()
    println("Error: You have not defined any input data to run the pipeline!")
    println("use --data <path(s)_to_data>")
    }

// check if the input is defined in a csv and parse the csv content to create the required channels
else if (params.data =~ /\.csv$/) {

    Channel
        .fromPath(params.data)
        .splitCsv(header:true)
        .map{ row-> tuple(file(row.receptor_path), file(row.ligand_path_sdf), file(row.ligand_path_mol2), file(row.reference_path), row.complex_name, row.BS ) }
        .branch{
            with_complex_name: it[4] != '' && it[4] != '-'
            other: true
        }
        .set{ all_input }

    // define a complex name where it is not given
    all_input.other.map{ [ it[0], it[1], it[2], it[3], it[0].simpleName + '__' + it[1].simpleName , it[5] ] }
                   .concat(all_input.with_complex_name)
                   .set{ all_input_defined }

    // map pdb, sdf, and mol2 files
    all_input_defined.map{ [it[0]] }.flatten().set{ pdb_files }
    all_input_defined.map{ [it[1]] }.flatten().set{ ref_sdf_files }
    all_input_defined.map{ [it[2]] }.filter{ it[].simpleName == /-/ }.ifEmpty( [] ).set{ mol_files }

    // ref files; if none is given assign the input pdb as reference
    all_input_defined.branch{
                        with_ref_file: it[3] != '' && it[3] != '-'
                        other: true
                        }
                        .set{ ref_definition }

    ref_definition.other.map{ [it[4], it[0].simpleName, it[1].simpleName, it[0], it[0], it[1]]}         // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
                        .concat(ref_definition.with_ref_file.map{ [it[4], it[0].simpleName, it[1].simpleName, it[3], it[0], it[1]]})
                        .set{ scoring_ref }

    // identifiers for receptor, ligand, complex names for each protein-ligand complex
    if (params.naming == "default") {
        //all_input_defined.map{ [it[0].simpleName, it[1].simpleName.split("__")[1].split("_")[0], it[4] ] }.set{ identifiers }
        all_input_defined.map{ [it[0].simpleName, it[1].simpleName.split("__")[1], it[4] ] }.set{ identifiers }
    }
    else {
        all_input_defined.map{ [it[0].simpleName, it[1].simpleName, it[4] ] }.set{ identifiers }
    }

    // get complexes with defined binding site
    all_input_defined.branch{
                        with_bs: it[5] != '' && it[5] != '-'
                        other: true
                        }
                        .set{ binding_sites }
    params.BS = true
}


// channel definitions if input files are in one directory
else if (params.data.split(',').size() == 1 ) {

    Channel
        .fromPath("${params.data}/*.sdf")
        .set { ref_sdf_files }

    Channel
        .fromPath("${params.data}/*.mol2")
        .set { mol_files }

    Channel
        .fromPath("${params.data}/*.pdb")
        .set { pdb_files }

    Channel
        .fromPath("${params.ref_files}/*.{pdb,cif}")
        .set { reference_files }

    if (params.naming == "default") {
        ref_sdf_files.map { [it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName] }
                     .set { identifiers }
    }
    else {
        println("Warning! Use the default naming scheme otherwise the receptor and the ligand might not be combined correctly!")
    }

    identifiers.combine(reference_files.map{ [it.simpleName, it] }, by: 0)            // receptor, ligand, complex, ref_receptor_file
               .combine(pdb_files.map{ [it.simpleName, it] }, by: 0)                // receptor, ligand, complex, ref_receptor_file, model_receptor_file
               .combine(ref_sdf_files.map{ [it , it.simpleName.split("__")[1]] }, by: 1)      // ligand, receptor, complex, ref_receptor_file, model_receptor_file, ref_ligand_file
               .map{ [it[2], it[1], it[0], it[3], it[4], it[5]] }.view()                   // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
               .set{ scoring_ref }

    params.BS = false
}

// channel definitions if input files are in two directories
else if (params.data.split(',').size() == 2 ) {

    Channel
        .fromPath("${params.data}".split(',')[1] + "/*.sdf")
        .set { ref_sdf_files }

    Channel
        .fromPath("${params.data}".split(',')[1] + "/*.mol2")
        .set { mol_files }

    Channel
        .fromPath("${params.data}".split(',')[0] + "/*.pdb")
        .set { pdb_files }

    Channel
        .fromPath("${params.ref_files}".split(',')[0] + "/*.{pdb,cif}")
        .set { reference_files }

    if (params.naming == "default") {
        ref_sdf_files.map { [it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName] }
                     .set { identifiers }
    }
    else {
        println("Warning! Use the default naming scheme otherwise the receptor and the ligand might not be combined correctly!")
    }

    identifiers.combine(reference_files.map{ [it.simpleName, it] }, by: 0)            // receptor, ligand, complex, ref_receptor_file
               .combine(pdb_files.map{ [it.simpleName, it] }, by: 0)                // receptor, ligand, complex, ref_receptor_file, model_receptor_file
               .combine(ref_sdf_files.map{ [it , it.simpleName.split("__")[1]] }, by: 1)      // ligand, receptor, complex, ref_receptor_file, model_receptor_file, ref_ligand_file
               .map{ [it[2], it[1], it[0], it[3], it[4], it[5]] }.view()                   // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
               .set{ scoring_ref }

    params.BS = false
}

// add Diffdock dependencies
Channel
    .fromPath("${params.diffdock_tool}/*", type: 'any')
    .set { diffd_tool }


/*
* include the modules
*/

include { ligand_preprocessing } from "./modules/ligand_preprocessing"
include { add_Hs_to_receptor } from "./modules/add_Hs_to_receptor"
include { p2rank } from "./modules/p2rank"
include { calculate_boxSize } from "./modules/calculate_boxSize"
include { docking_boxes_predicted_pockets; docking_box_defined_BS } from "./modules/docking_box"
include { diffdock; diffdock_single } from "./modules/diffdock"
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

    // ligand preprocessing
    sdf_for_docking = ligand_preprocessing(ref_sdf_files.unique().collect(), mol_files.unique().collect().ifEmpty([]))

    sdf_for_docking.sdf_files.flatten()
                             .map{ [it.simpleName.toString().split("_preped")[0], it] }
                             .combine(identifiers.map{ it[1] }, by: 0)
                             .unique()
                             .set { ligand_tuple }
    ligand_tuple.map{ it[1] }.set { ligand_only }

    // add Hs to receptor
    if (params.receptor_Hs == "no") {
        pdb_Hs = add_Hs_to_receptor(pdb_files.unique().map{ [it.simpleName, it] } )
    }
    else {
        pdb_files.unique().map { [it.simpleName, it] }.set{ pdb_Hs }
    }


    /*
    * binding pocket prediction
    */

    // predict binding pockets
    binding_pockets = p2rank(pdb_Hs)

    // define box parameters for vina-like tools
    if (params.tools =~ /vina/ || params.tools =~ /smina/ || params.tools =~ /gnina/) {
        box_size = calculate_boxSize(ligand_tuple)

        if (params.BS) {
            // boxes without defined BS
            binding_sites.other.map{ [it[4]] }                               // complex
                               .combine(identifiers.map{ [it[2], it[0], it[1]] }, by: 0).map{ [it[1], it[2], it[0]] }     // receptor, ligand, complex
                               .combine(pdb_Hs, by: 0)                      // receptor, ligand, complex, pdb_Hs
                               .combine(binding_pockets.pockets, by: 0)     // receptor, ligand, complex, pdb_Hs, binding_pockets
                               .combine(box_size.size, by: 1)               // ligand, receptor, complex, pdb_Hs, binding_pockets, box_size
                               .set{ input_dockingBox_no_BS }

            boxes_no_BS = docking_boxes_predicted_pockets(input_dockingBox_no_BS)

            // boxes with defined BS
            binding_sites.with_bs.map{ [it[4], it[5]] }                     // complex, BS
                                 .combine(identifiers.map{ [it[2], it[0], it[1]] }, by: 0).map{ [it[2], it[3], it[0], it[1]] }     // receptor, ligand, complex, BS
                                 .combine(pdb_Hs, by: 0)                    // receptor, ligand, complex, BS, pdb_Hs
                                 .combine(box_size.size, by: 1)             // ligand, receptor, complex, BS, pdb_Hs, box_size
                                 .map{ [it[0], it[1], it[2], it[4], it[3], it[5]] }      // ligand, receptor, complex, pdb_Hs, BS, box_size
                                 .set{ input_dockingBox_with_BS }

            boxes_BS = docking_box_defined_BS(input_dockingBox_with_BS)

            // concatenate boxes_no_BS and boxes_BS channels
            boxes = boxes_no_BS.concat(boxes_BS)
        }
        else {
            identifiers.combine(pdb_Hs, by: 0)
                       .combine(binding_pockets.pockets, by: 0)
                       .combine(box_size.size, by: 1)
                       .set{ input_dockingBox }
            boxes = docking_boxes_predicted_pockets(input_dockingBox)
        }
    }


    /*
    * docking using Diffdock
    */

    if (params.tools =~ /diffdock/) {

        if (params.diffdock_mode == "batch") {

            identifiers.combine(pdb_Hs, by: 0)
                       .combine(ligand_tuple.map{ [it[1], it[0]]}, by: 1)
                       .map{ [it[2], it[3].name.toString(), it[4].name.toString()] }
                       .collectFile() {item ->
                            [ "protein_ligand.csv", item[0] + "," + item[1] + "," + item[2] + ",\n"]
                       }
                       .set{ diffd_csv }

            diffdock_predictions = diffdock(diffd_csv, pdb_Hs.flatten().filter{it =~ /\//}.collect(), sdf_for_docking.sdf_files.collect(), diffd_tool.collect())
        }

        else if (params.diffdock_mode == "single") {

            identifiers.combine(pdb_Hs, by: 0)
                       .combine(ligand_tuple.map{ [ it[1], it[0] ] }, by: 1)
                       .set{ input_diffd_single }

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

        identifiers.combine(preped_receptors, by: 0)                                // receptor, ligand, complex, preped_receptor
                   .combine(preped_ligands.map{ [it[1], it[0]] }, by: 1)            // ligand, receptor, complex, preped_receptor, preped_ligand
                   .map { [ it[2], it[0], it[1], it[3], it[4] ] }                   // complex, ligand, receptor, preped_receptor, preped_ligand
                   .combine(boxes.flatten().map{file -> tuple(file.simpleName.toString().split("_pocket")[0], file)}, by: 0)    // complex, ligand, receptor, preped_receptor, preped_ligand, box_file
                   .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }            // complex, ligand, receptor, pocket, preped_receptor, preped_ligand, box_file
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

    // scoring the modelled receptors

    if (params.alphafold == "yes") {
        modelled_receptors_scores = ost_scoring_modelReceptors(scoring_ref.map{ [it[1], it[3], it[4]] }.unique())
        modelled_receptors_scores_summary = combine_modelReceptors_scores(modelled_receptors_scores.toList().flatten().filter{ it =~ /json/ }.collect())
    }


    // tankbind

    if (params.tools =~ /tankbind/) {
        scoring_ref.combine(tankbind_out.sdfs, by: 0)
                   .set { tb_scoring_input }
        tb_scores = tb_ost(tb_scoring_input, Channel.value( 'tankbind' ))
        tb_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ tb_scores_for_summary }
    }


    // diffdock
    if (params.tools =~ /diffdock/) {
        scoring_ref.combine(diffdock_predictions.predictions.flatten().filter{ it =~ /\// }.map{[it.toString().split('/')[-2], it]}.groupTuple(), by: 0)
                   .set { dd_scoring_input }
        dd_scores = dd_ost(dd_scoring_input, Channel.value( 'diffdock' ))
        dd_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ dd_scores_for_summary }
    }


    // vina
    if (params.tools =~ /vina/) {
        scoring_ref.combine(vina_sdf.map{[ it[0], it[3]]}, by: 0)
                   .set { vina_scoring_input }
        vina_scores = vina_ost(vina_scoring_input, Channel.value( 'vina' ))
        vina_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ vina_scores_for_summary }
    }


    // smina
    if (params.tools =~ /smina/) {
        scoring_ref.combine(smina_out.smina_sdf.map{[ it[0], it[3]]}, by: 0)
                   .set { smina_scoring_input }
        smina_scores = smina_ost(smina_scoring_input, Channel.value( 'smina' ))
        smina_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ smina_scores_for_summary }
    }


    // gnina
    if (params.tools =~ /gnina/) {
        scoring_ref.combine(gnina_out.gnina_sdf.map{[ it[0], it[3]]}, by: 0)
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
