#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Define the pipeline parameters
 */

log.info """
PickyBinder  ~  version ${workflow.manifest.version}
=============================================

Input data             : ${params.data}
Reference files        : ${params.ref_files}

Tools                  : ${params.tools}
Diffdock mode          : ${params.diffdock_mode}

Ligand scoring         : ${params.scoring_ligands}
Receptor scoring       : ${params.scoring_receptors}

Hydrated receptors     : ${params.receptor_Hs}
Autobox add            : ${params.autobox_add}

=============================================

"""


/*
* define the input channels
*/

// check if there is a path to the input data
if ( params.data == "" ) {
    println()
    println()
    println("Error: You have not defined any input data to run the pipeline!")
    println("use --data <path(s)_to_data>")
    }

// check if the input is defined in a csv and parse the csv content to create the required channels
else if ( params.data =~ /\.csv$/ ) {

    Channel
        .fromPath( params.data )
        .splitCsv( header:true )
        .map{ row -> tuple( file(row.receptor), row.ligand, file(row.ligand_mol2),
                            file(row.reference), row.complex_name, row.BS, row.alphafold ) }
        .branch{
            ligand_file: it[1].endsWith(".sdf") || it[1].endsWith(".smi") || it[1].endsWith(".mol2")
            other: true
        }
        .set{ ligand_type }

    ligand_type.other.collectFile() { item -> [ "${item[4]}.smi", item[1] + '\n' ] }
                     .set{ smi_files }

    ligand_type.other.map{ [ it[4], it[0], it[1], it[2], it[3], it[5], it[6] ] } // complex, pdb_file, smiles_string, mol2_file, ref_file, BS, alphafold
                     .combine( smi_files.map{ [ it.simpleName, it ] }, by: 0 ) // complex [0], pdb_file [1], smiles_string [2], mol2_file [3], ref_file [4], BS [5], alphafold [6], smi_file [7]
                     .map{ [ it[1], it[7], it[3], it[4], it[0], it[5], it[6] ] } // pdb_file, smi_file, mol2_file, ref_file, complex, BS, alphafold
                     .concat( ligand_type.ligand_file.map{ [ it[0], file( it[1] ), it[2], it[3], it[4], it[5], it[6] ] } )
                     .branch{
                        with_complex_name: it[4] != '' && it[4] != '-'
                        other: true
                     }
                    .set{ all_input }

    // define a complex name where it is not given
    all_input.other.map{ [ it[0], it[1], it[2], it[3], it[0].simpleName + '__' + it[1].simpleName , it[5],  it[6] ] }
                   .concat( all_input.with_complex_name )
                   .set{ all_input_defined }

    // map pdb, sdf, and mol2 files
    all_input_defined.map{ [ it[0], it[6] ] }.set{ pdb_files }
    all_input_defined.map{ [ it[1] ] }.flatten().set{ ref_sdf_files }
    all_input_defined.map{ [ it[2] ] }.filter{ !(it.simpleName == /-/) }.ifEmpty( [] ).set{ mol_files }

    // ref files; if none is given assign the input pdb as reference
    all_input_defined.branch{
                        with_ref_file: it[3].simpleName != '' && it[3].simpleName != '-'
                        other: true
                        }
                        .set{ ref_definition }

    ref_definition.other.map{ [ it[4], it[0].simpleName, it[1].simpleName, it[0], it[0], it[1] ] }   // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
                        .concat( ref_definition.with_ref_file.map{ [ it[4], it[0].simpleName, it[1].simpleName, it[3], it[0], it[1] ] } )
                        .set{ scoring_ref }

    // identifiers for receptor, ligand, complex names for each protein-ligand complex
    all_input_defined.map{ [ it[0].simpleName, it[1].simpleName, it[4] ] }
                     .set{ identifiers }

    // get complexes with defined binding site
    all_input_defined.branch{
                        with_bs: it[5] != '' && it[5] != '-'
                        other: true
                        }
                        .set{ binding_sites }
    params.BS = true
}


// channel definitions if input files are in one directory
else if ( params.data.split(',').size() == 1 ) {

    Channel
        .fromPath( "${params.data}/*.sdf" )
        .set { ref_sdf_files }

    Channel
        .fromPath( "${params.data}/*.mol2" )
        .set { mol_files }

    if ( params.alphafold == "yes" )
        Channel
            .fromPath( "${params.data}/*.pdb" )
            .map{ [ it, "yes" ] }
            .set { pdb_files }
    else {
        Channel
            .fromPath( "${params.data}/*.pdb" )
            .map{ [ it, "" ] }
            .set { pdb_files }
    }

    Channel
        .fromPath( "${params.ref_files}/*.{pdb,cif}" )
        .set { reference_files }

    ref_sdf_files.map { [ it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName ] }
                 .set { identifiers }

    identifiers.combine( reference_files.map{ [ it.simpleName, it ] }, by: 0 )          // receptor, ligand, complex, ref_receptor_file
               .combine( pdb_files.map{ [it[0].simpleName, it[0] ] }, by: 0 )                // receptor, ligand, complex, ref_receptor_file, model_receptor_file
               .combine( ref_sdf_files.map{ [ it , it.simpleName.split("__")[1] ] }, by: 1 )      // ligand, receptor, complex, ref_receptor_file, model_receptor_file, ref_ligand_file
               .map{ [ it[2], it[1], it[0], it[3], it[4], it[5] ] }                   // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
               .set{ scoring_ref }

    params.BS = false
}

// channel definitions if input files are in two directories
else if ( params.data.split(',').size() == 2 ) {

    Channel
        .fromPath( "${params.data}".split(',')[1] + "/*.sdf" )
        .set { ref_sdf_files }

    Channel
        .fromPath( "${params.data}".split(',')[1] + "/*.mol2" )
        .set { mol_files }

    if ( params.alphafold == "yes" )
        Channel
            .fromPath( "${params.data}".split(',')[0] + "/*.pdb" )
            .map{ [ it, "yes" ] }
            .set { pdb_files }
    else {
        Channel
            .fromPath( "${params.data}".split(',')[0] + "/*.pdb" )
            .map{ [ it, "" ] }
            .set { pdb_files }
    }

    Channel
        .fromPath( "${params.ref_files}".split(',')[0] + "/*.{pdb,cif}" )
        .set { reference_files }

    ref_sdf_files.map { [ it.simpleName.split("__")[0], it.simpleName.split("__")[1], it.simpleName ] }
                 .set { identifiers }

    identifiers.combine( reference_files.map{ [ it.simpleName, it ] }, by: 0 )          // receptor, ligand, complex, ref_receptor_file
               .combine( pdb_files.map{ [ it[0].simpleName, it[0] ] }, by: 0 )                // receptor, ligand, complex, ref_receptor_file, model_receptor_file
               .combine( ref_sdf_files.map{ [ it , it.simpleName.split("__")[1] ] }, by: 1 )      // ligand, receptor, complex, ref_receptor_file, model_receptor_file, ref_ligand_file
               .map{ [ it[2], it[1], it[0], it[3], it[4], it[5] ] }                   // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
               .set{ scoring_ref }

    params.BS = false
}

// add Diffdock and EDMDock dependencies
Channel
    .fromPath( "${params.diffdock_tool}/*", type: 'any' )
    .set { diffd_tool }

Channel
    .fromPath( "${params.edmdock_tool}/*", type: 'any' )
    .set { edmdock_tool }


/*
* include the modules
*/

include { ligand_preprocessing_single; ligand_preprocessing_log } from "./modules/ligand_preprocessing"
include { add_Hs_to_receptor; fix_pdb } from "./modules/add_Hs_to_receptor"
include { p2rank } from "./modules/p2rank"
include { calculate_boxSize } from "./modules/calculate_boxSize"
include { docking_boxes_predicted_pockets; docking_box_defined_BS } from "./modules/docking_box"
include { diffdock; diffdock_single } from "./modules/diffdock"
include { vina_prepare_receptor; vina_prepare_ligand; vina; pdbtqToSdf as vina_pdbtqToSdf; pdbtqToSdf as smina_pdbtqToSdf; pdbtqToSdf as gnina_pdbtqToSdf } from "./modules/vina"
include { gnina_sdf } from "./modules/gnina"
include { smina_sdf } from "./modules/smina"
include { tankbind } from "./modules/tankbind"
include { edmdock; edmdock_single } from "./modules/edmdock"
include { pdb_to_sdf_batch as edmdock_pdb_to_sdf_batch; pdb_to_sdf_single as edmdock_pdb_to_sdf_single} from "./modules/scoring"
include { ost_scoring as tb_ost; ost_scoring as dd_ost; ost_scoring as vina_ost; ost_scoring as smina_ost; ost_scoring as gnina_ost; ost_scoring as edm_ost; ost_score_summary } from "./modules/scoring"
include { ost_scoring_receptors; combine_receptors_scores } from "./modules/scoring"
include { combine_all_scores } from "./modules/all_scores_summary"
include { catch_ignored_tasks; catch_diffdock_problems; no_box_size; error_and_problems_summary } from "./modules/catch_failed_tasks"


/*
* main workflow
*/

workflow {

    /*
    * receptor preprocessing
    */

    // add Hs to receptor
    if ( params.receptor_Hs == "no" ) {
        pdb_Hs = add_Hs_to_receptor( pdb_files.unique().map{ [ it[0].simpleName, it[0], it[1] ] } )
    }
    else {
        pdb_files.unique().map{ [ it[0].simpleName, it[0], it[1] ] }
                          .set{ pdb_Hs }
    }

    /*
    * binding pocket prediction
    */

    // predict binding pockets
    if ( params.tools =~ /vina/ || params.tools =~ /smina/ || params.tools =~ /gnina/ || params.tools =~ /tankbind/ || params.tools =~ /edmdock/ || params.p2rank_only == "yes" ) {
        binding_pockets = p2rank( pdb_Hs )
        params.P2RANK = "done"
    }
    else {
        params.P2RANK = "no"
    }

    if ( params.p2rank_only != "yes" )  {

        /*
        * ligand preprocessing
        */

        // ligand preprocessing
        sdf_for_docking = ligand_preprocessing_single( all_input_defined.map{ [ it[1].simpleName, it[1], it[2] ] } )

        ligand_prep_log = ligand_preprocessing_log( sdf_for_docking.ligand_prep_log.collect() )

        sdf_for_docking.sdf_files.flatten()
                                 .map{ [ it.simpleName.toString().split("_preped")[0], it ] }
                                 .combine( identifiers.map{ it[1] }, by: 0 )
                                 .unique()
                                 .set { ligand_tuple }
        ligand_tuple.map{ it[1] }.set { ligand_only }


        // define box parameters for vina-like tools
        wp_coordinates = Channel.empty()
        if ( params.tools =~ /vina/ || params.tools =~ /smina/ || params.tools =~ /gnina/ || params.tools =~ /edmdock/ ) {
            box_size = calculate_boxSize( ligand_tuple )
            box_size_failed = no_box_size( box_size.log.collect() )

            if (params.BS) {
                // boxes without defined BS
                binding_sites.other.map{ [ it[4] ] }                                     // complex
                                   .combine( identifiers.map{ [ it[2], it[0], it[1] ] }, by: 0 ).map{ [ it[1], it[2], it[0] ] }  // receptor, ligand, complex
                                   .combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )     // receptor, ligand, complex, pdb_Hs
                                   .combine( binding_pockets.pockets, by: 0 )            // receptor, ligand, complex, pdb_Hs, binding_pockets
                                   .combine( box_size.size, by: 1 )                      // ligand, receptor, complex, pdb_Hs, binding_pockets, box_size
                                   .set{ input_dockingBox_no_BS }

                boxes_no_BS = docking_boxes_predicted_pockets( input_dockingBox_no_BS )

                // boxes with defined BS
                binding_sites.with_bs.map{ [ it[4], it[5] ] }                              // complex, BS
                                     .combine( identifiers.map{ [ it[2], it[0], it[1] ] }, by: 0 ).map{ [ it[2], it[3], it[0], it[1] ] }  // receptor, ligand, complex, BS
                                     .combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )     // receptor, ligand, complex, BS, pdb_Hs
                                     .combine( box_size.size, by: 1 )                      // ligand, receptor, complex, BS, pdb_Hs, box_size
                                     .map{ [ it[0], it[1], it[2], it[4], it[3], it[5] ] }  // ligand, receptor, complex, pdb_Hs, BS, box_size
                                     .set{ input_dockingBox_with_BS }

                boxes_BS = docking_box_defined_BS( input_dockingBox_with_BS )

                // concatenate boxes_no_BS and boxes_BS channels
                boxes = boxes_no_BS.box_per_pocket.concat( boxes_BS.boxes_bs_wp )
                wp_coordinates = boxes_no_BS.center_coordinates.concat( boxes_BS.center_coordinates )
            }
            else {
                identifiers.combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )
                           .combine( binding_pockets.pockets, by: 0 )
                           .combine( box_size.size, by: 1 )
                           .set{ input_dockingBox }
                boxes_all = docking_boxes_predicted_pockets( input_dockingBox )
                boxes_all.box_per_pocket.set{ boxes }
                boxes_all.center_coordinates.set{ wp_coordinates }
            }
        }
        else {
            box_size_failed = Channel.empty()
        }


        /*
        * docking using Diffdock
        */

        if ( params.tools =~ /diffdock/ ) {

            if ( params.diffdock_mode == "batch" ) {

                identifiers.combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )
                           .combine( ligand_tuple.map{ [ it[1], it[0] ]}, by: 1 )
                           .map{ [ it[2], it[3].name.toString(), it[4].name.toString() ] }

                           .collectFile() {item ->
                                [ "protein_ligand.csv", item[0] + "," + item[1] + "," + item[2] + ",\n"]
                           }
                           .set{ diffd_csv }

                diffdock_predictions = diffdock( diffd_csv, pdb_Hs.map{ [ it[0], it[1] ] }.flatten().filter{ it =~ /\// }.collect(), sdf_for_docking.sdf_files.collect(), diffd_tool.collect() )
                dd_problems = catch_diffdock_problems( diffdock_predictions.log )

                Channel.fromPath( "$launchDir/predictions/diffdock/diffdock_predictions/*/*confidence*" ).ifEmpty( [] )
                       .concat( diffdock_predictions.predictions.flatten().ifEmpty( [] ))
                       .filter{ it =~ /confidence/ }
                       .unique()
                       .set{ all_dd_predictions }

                all_dd_predictions.collectFile() {item -> [ "dd_file_names.txt", item.name + "\n" ] }
                                  .set{ for_dd_tool_scores }
            }
            else if ( params.diffdock_mode == "single" ) {

                identifiers.combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )
                           .combine( ligand_tuple.map{ [ it[1], it[0] ] }, by: 1 )
                           .set{ input_diffd_single }

                diffdock_predictions = diffdock_single( input_diffd_single, diffd_tool.collect() )

                diffdock_predictions.predictions.flatten()
                                    .filter{ it =~ /confidence/ }
                                    .set{ all_dd_predictions }

                all_dd_predictions.collectFile() {item -> [ "dd_file_names.txt", item.name + "\n" ] }
                                  .set{ for_dd_tool_scores }

                dd_problems = Channel.empty()
            }
        }
        else {
            dd_scores_for_summary = Channel.empty()
            for_dd_tool_scores = Channel.empty()
            dd_problems = Channel.empty()
        }


        /*
        * docking using Vina
        */

        if ( params.tools =~ /vina/ ) {
            preped_ligands = vina_prepare_ligand( ligand_tuple )
            preped_receptors = vina_prepare_receptor( pdb_Hs.map{ [ it[0], it[1] ] } )

            identifiers.combine( preped_receptors, by: 0 )                        // receptor, ligand, complex, preped_receptor
                       .combine( preped_ligands.map{ [ it[1], it[0] ] }, by: 1 )    // ligand, receptor, complex, preped_receptor, preped_ligand
                       .map { [ it[2], it[0], it[1], it[3], it[4] ] }           // complex, ligand, receptor, preped_receptor, preped_ligand
                       .combine( boxes.transpose(), by: 0 )                       // complex, ligand, receptor, preped_receptor, preped_ligand, box_file
                       .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }     // complex, ligand, receptor, pocket, preped_receptor, preped_ligand, box_file
                       .set { vina_input }

            vina_out = vina( vina_input )
            vina_sdf = vina_pdbtqToSdf( vina_out.vina_result, Channel.value( 'vina' ) )
            vina_sdf.map{ it[3] }.collect().set{ for_vina_tool_scores }
        }
        else {
            vina_scores_for_summary = Channel.empty()
            for_vina_tool_scores = Channel.empty()
        }


        /*
        * docking using smina
        */

        if ( params.tools =~ /smina/ || params.tools =~ /gnina/ ) {
            identifiers.combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )
                       .combine( ligand_tuple.map{ [ it[1], it[0] ] }, by: 1 )
                       .map { [ it[2], it[0], it[1], it[3], it[4] ] }
                       .combine( boxes.transpose(), by: 0 )
                       .map { [ it[0], it[1], it[2], it[5].simpleName.toString().split("_")[-1], it[3], it[4], it[5] ] }
                       .set { smi_gni_input }
        }

        if ( params.tools =~ /smina/ ) {
            smina_out = smina_sdf( smi_gni_input )
            smina_out.smina_sdf.map{ it[3] }.collect().set{ for_smina_tool_scores }
        }
        else {
            smina_scores_for_summary = Channel.empty()
            for_smina_tool_scores = Channel.empty()
        }


        /*
        * docking using gnina
        */

        if ( params.tools =~ /gnina/ ) {
            gnina_out = gnina_sdf( smi_gni_input )
            gnina_out.gnina_sdf.map{ it[3] }.collect().set{ for_gnina_tool_scores }
        }
        else {
            gnina_scores_for_summary = Channel.empty()
            for_gnina_tool_scores = Channel.empty()
        }


        /*
        * docking using tankbind
        */

        if ( params.tools =~ /tankbind/ ) {
            identifiers.combine( pdb_Hs.map{ [ it[0], it[1] ] }, by: 0 )
                       .combine( binding_pockets.pockets, by: 0 )
                       .combine( ligand_only.map{ [ it, it.simpleName.toString().split("_preped")[0] ] }, by: 1 )
                       .set{ tankbind_input }

            tankbind_out = tankbind( tankbind_input )
            tankbind_out.affinities.map{ it[2] }.collect().set{ for_tb_tool_scores }
        }
        else {
            tb_scores_for_summary = Channel.empty()
            for_tb_tool_scores = Channel.empty()
        }


        /*
        * docking using EDMDock
        */

        if ( params.tools =~ /edmdock/ ) {

            pdbfixer_fixed_pdbs = fix_pdb( pdb_Hs.map{ [ it[0], it[1] ] } )

            identifiers.combine( pdbfixer_fixed_pdbs, by: 0 )
                       .combine( ligand_tuple.map{ [ it[1], it[0] ] }, by: 1 )
                       .map { [ it[2], it[0], it[1], it[3], it[4] ] }
                       .combine( boxes.transpose(), by: 0 )   // complex, ligand, receptor, preped_receptor, preped_ligand, box_file
                       .map { [ it[0], it[5].simpleName.toString().split("_")[-1], it[2], it[1], it[3], it[4], it[5] ] }  // complex, pocket, receptor, ligand, preped_receptor, preped_ligand, box_file
                       .set { edm_dock_input }

            // single sample version
            edmdock_out = edmdock_single( edm_dock_input, edmdock_tool.collect() )
            edmdock_sdfs = edmdock_pdb_to_sdf_single( edmdock_out.map{ [ it[0], it[1], it[3] ] }, "predictions/edmdock/sdf_files" )
            for_edmdock_tool_scores = edmdock_sdfs.collect()
        }
        else {
            edmdock_scores_for_summary = Channel.empty()
            for_edmdock_tool_scores = Channel.empty()
        }


        /*
        * ost scoring
        */

        //get coordinates used for docking box definition

        if ( params.BS ) {
            // coordinates whole protein box
            wp_coordinates.map{ complex, receptor, csv -> [ complex, receptor, csv.splitCsv(header:true, strip:true) ] }
                          .map{ complex, receptor, row -> [ receptor, "pocketProteinCenter", complex, row.center_x, row.center_y, row.center_z ] }
                          .transpose()
                          .set{ receptor_coordinates }        // receptor, pocket, complex, x, y, z

            // defined BS coordinates
            binding_sites.with_bs.map{ [ it[0].simpleName, "pocketBS", it[4], it[5].split("_")[0], it[5].split("_")[1], it[5].split("_")[2] ] }
                                 .set{ bs_coordinates }     //  receptor, pocket, complex, x, y, z

            bs_coordinates.concat( receptor_coordinates ).set{ defined_coordinates }      //  receptor, pocket, complex, x, y, z
        }
        else {
            // coordinates whole protein box
            wp_coordinates.map{ complex, receptor, csv -> [ complex, receptor, csv.splitCsv(header:true, strip:true) ] }
                          .map{ complex, receptor, row -> [ receptor, "pocketProteinCenter", complex, row.center_x, row.center_y, row.center_z ] }
                          .transpose()
                          .set{ defined_coordinates }        // receptor, pocket, complex, x, y, z
        }

        // coordinates from p2rank predictions
        if ( params.tools =~ /vina/ || params.tools =~ /smina/ || params.tools =~ /gnina/ || params.tools =~ /tankbind/ || params.tools =~ /edmdock/ ) {
            binding_pockets.pockets.map{ receptor, csv -> [ receptor, csv.splitCsv(header:true, strip:true) ] }
                                   .branch{
                                        p2rank_ok: it[1] != []
                                        other: true
                                   }
                                   .set{ p2rank_success }
            p2rank_success.p2rank_ok.map{ receptor, row -> [ receptor, row.name, row.center_x, row.center_y, row.center_z ] }
                                    .transpose()
                                    .set{ predicted_coordinates }  // receptor, pocket, x, y, z

            predicted_coordinates.map { it -> "${it[0]},${it[1]},${it[2]},${it[3]},${it[4]}" }
                                 .collectFile( name: 'p2rank_summary.csv', newLine: true )
                                 .set { p2rank_summary }

            p2rank_success.other.collectFile( storeDir: 'preprocessing/p2rank' ) {item -> [ "p2rank_no_pockets_found.csv", item[0] + "\n"] }
                                .set{ p2rank_no_pocket }
        }
        else {
            p2rank_no_pocket = Channel.empty()
        }

        // scoring the modelled receptors
        if ( params.scoring_receptors == "yes" ) {
            receptors_scores = ost_scoring_receptors( scoring_ref.map{ [ it[1], it[3], it[4] ] }.unique() )
            receptors_scores_summary = combine_receptors_scores( receptors_scores.toList().flatten().filter{ it =~ /json/ }.collect() )
        }

        if ( params.scoring_ligands == "yes" ) {

            // tankbind

            if ( params.tools =~ /tankbind/ ) {
                scoring_ref.combine( tankbind_out.sdfs, by: 0 )
                           .set{ tb_scoring_input }
                tb_scores = tb_ost( tb_scoring_input, Channel.value( 'tankbind' ) )

                tb_scores.summary.map{ complex, receptor, csv -> [ complex, receptor, csv.splitCsv(header:true, strip:true) ] }
                                 .map{ complex, receptor, row -> [ complex, receptor, row.Tool, row.Complex, row.Pocket, row.Rank, row.'lDDT-PLI', row.'lDDT-LP', row.BiSyRMSD, row.Reference_Ligand ] }
                                 .transpose()
                                 .collectFile() { item -> [ "${item[1]}____${item[0]}_${item[4]}_tankbind_score_summary.csv", item[2] + "," + item[3] + "," + item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] + "," + item[9] ] }
                                 .map{ [ it.simpleName.split('____')[0], it.simpleName.split('_')[-4], it.simpleName.split('_pocket')[0].split('____')[1], it ] }
                                 .set{ tb_scores_for_coordinates }    // receptor, pocket, complex, tb_score_csv

                tankbind_out.affinities.map{ complex, receptor, csv -> [ complex, receptor, csv.splitCsv(strip:true, skip: 1) ] }
                                       .transpose()
                                       .map { [ it[1], it[2][2], it[0], it[2][3].split('"')[1], it[2][4], it[2][5].split('"')[0] ] }    // receptor, pocket, complex, x, y, z
                                       .combine( tb_scores_for_coordinates, by: [0, 1, 2] )
                                       .map{ receptor, pocket, complex, x, y, z, csv -> [ receptor, pocket, complex, x, y, z, csv.splitCsv(strip:true)] }
                                       .transpose()
                                       .map { [ it[0], it[1], it[2], it[6][0], it[6][1], it[6][2], it[6][3], it[6][4], it[6][5], it[6][6], it[6][7], it[3], it[4], it[5] ] }
                                       .collectFile() { item -> [ "${item[2]}_${item[1]}_tankbind_score_summary.csv", item[3] + "," + item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] + "," + item[9] + "," + item[10] + "," + item[11] + "," + item[12] + "," + item[13] ] }
                                       .collect()
                                       .set{ tb_scores_for_summary }
            }


            // diffdock
            if ( params.tools =~ /diffdock/ ) {
                scoring_ref.combine( all_dd_predictions.map{ [ it.toString().split('/')[-2], it ] }.groupTuple(), by: 0 )
                           .set { dd_scoring_input }
                dd_scores = dd_ost( dd_scoring_input, Channel.value( 'diffdock' ) )
                dd_scores.summary.toList().flatten().filter{ it =~ /\.csv/ }.collect().set{ dd_scores_for_summary }
            }


            // vina
            if ( params.tools =~ /vina/ ) {
                scoring_ref.combine( vina_sdf.map{ [ it[0], it[3] ] }, by: 0 )
                           .set { vina_scoring_input }
                vina_scores = vina_ost( vina_scoring_input, Channel.value( 'vina' ) )

                vina_scores.summary.map{ [ it[1], "pocket" + it[2].simpleName.toString().split("_pocket")[1].split("_")[0], it[0], it[2] ] }   // receptor, pocket, complex, vina_score_csv
                                   .set{ vina_scores_for_coordinates }

                vina_scores_for_coordinates.combine( predicted_coordinates, by: [0, 1] )        // receptor, pocket, complex, vina_score_csv, x, y, z
                                           .set{ vina_scores_for_coordinates_p2rank }

                vina_scores_for_coordinates.combine( defined_coordinates, by:  [0, 1, 2] )
                                           .set{ vina_scores_for_coordinates_defined }

                vina_scores_for_coordinates_p2rank.concat( vina_scores_for_coordinates_defined )
                            .map{ [ it[3].simpleName, it[0], it[1], it[2], it[3], it[4], it[5], it[6] ] }
                            .map{ file_name, receptor, pocket, complex, csv, x, y, z -> [ file_name, receptor, pocket, complex, csv.splitCsv(header:true, strip:true), x, y, z] }
                            .map{ file_name, receptor, pocket, complex, row, x, y, z -> [ file_name, receptor, pocket, complex, row.Tool, row.Complex, row.Pocket, row.Rank, row.'lDDT-PLI', row.'lDDT-LP', row.BiSyRMSD, row.Reference_Ligand, x, y, z] }
                            .transpose()
                            .collectFile() { item -> [ "${item[0]}.csv", item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] +  "," + item[9] + "," + item[10] + "," + item[11] + "," + item[12] + "," + item[13] + "," + item[14] + '\n'] }
                            .collect()
                            .set{ vina_scores_for_summary}
            }


            // smina
            if ( params.tools =~ /smina/ ) {
                scoring_ref.combine( smina_out.smina_sdf.map{ [ it[0], it[3] ] }, by: 0 )
                           .set { smina_scoring_input }
                smina_scores = smina_ost( smina_scoring_input, Channel.value( 'smina' ) )

                smina_scores.summary.map{ [ it[1], "pocket" + it[2].simpleName.toString().split("_pocket")[1].split("_")[0], it[0], it[2] ] }   // receptor, pocket, complex, vina_score_csv
                                    .set{ smina_scores_for_coordinates }

                smina_scores_for_coordinates.combine( predicted_coordinates, by: [0, 1] )        // receptor, pocket, complex, smina_score_csv, x, y, z
                                            .set{ smina_scores_for_coordinates_p2rank }

                smina_scores_for_coordinates.combine( defined_coordinates, by:  [0, 1, 2] )
                                            .set{ smina_scores_for_coordinates_defined }

                smina_scores_for_coordinates_p2rank.concat( smina_scores_for_coordinates_defined )
                         .map{ [ it[3].simpleName, it[0], it[1], it[2], it[3], it[4], it[5], it[6] ] }
                         .map{ file_name, receptor, pocket, complex, csv, x, y, z -> [ file_name, receptor, pocket, complex, csv.splitCsv(header:true, strip:true), x, y, z] }
                         .map{ file_name, receptor, pocket, complex, row, x, y, z -> [ file_name, receptor, pocket, complex, row.Tool, row.Complex, row.Pocket, row.Rank, row.'lDDT-PLI', row.'lDDT-LP', row.BiSyRMSD, row.Reference_Ligand, x, y, z] }
                         .transpose()
                         .collectFile() { item -> [ "${item[0]}.csv", item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] +  "," + item[9] + "," + item[10] + "," + item[11] + "," + item[12] + "," + item[13] + "," + item[14] + '\n'] }
                         .collect()
                         .set{ smina_scores_for_summary}
            }


            // gnina
            if ( params.tools =~ /gnina/ ) {
                scoring_ref.combine( gnina_out.gnina_sdf.map{ [ it[0], it[3] ] }, by: 0 )
                           .set { gnina_scoring_input }
                gnina_scores = gnina_ost( gnina_scoring_input, Channel.value( 'gnina' ) )

                gnina_scores.summary.map{ [ it[1], "pocket" + it[2].simpleName.toString().split("_pocket")[1].split("_")[0], it[0], it[2] ] }   // receptor, pocket, complex, vina_score_csv
                                    .set{ gnina_scores_for_coordinates }

                gnina_scores_for_coordinates.combine( predicted_coordinates, by: [0, 1] )        // receptor, pocket, complex, gnina_score_csv, x, y, z
                                            .set{ gnina_scores_for_coordinates_p2rank }

                gnina_scores_for_coordinates.combine( defined_coordinates, by:  [0, 1, 2] )
                                            .set{ gnina_scores_for_coordinates_defined }

                gnina_scores_for_coordinates_p2rank.concat( gnina_scores_for_coordinates_defined )
                         .map{ [ it[3].simpleName, it[0], it[1], it[2], it[3], it[4], it[5], it[6] ] }
                         .map{ file_name, receptor, pocket, complex, csv, x, y, z -> [ file_name, receptor, pocket, complex, csv.splitCsv(header:true, strip:true), x, y, z] }
                         .map{ file_name, receptor, pocket, complex, row, x, y, z -> [ file_name, receptor, pocket, complex, row.Tool, row.Complex, row.Pocket, row.Rank, row.'lDDT-PLI', row.'lDDT-LP', row.BiSyRMSD, row.Reference_Ligand, x, y, z] }
                         .transpose()
                         .collectFile() { item -> [ "${item[0]}.csv", item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] +  "," + item[9] + "," + item[10] + "," + item[11] + "," + item[12] + "," + item[13] + "," + item[14] + '\n'] }
                         .collect()
                         .set{ gnina_scores_for_summary}
            }


            // edmdock
            if ( params.tools =~ /edmdock/ ) {
                scoring_ref.combine( edmdock_sdfs.sdf_files.collect().flatten().map{ [ it.simpleName.split('_pocket')[0], it ] }.groupTuple(), by: 0 )
                           .set { edm_scoring_input }

                edm_scores = edm_ost( edm_scoring_input, Channel.value( 'edmdock' ) )

                edm_scores.summary.map{ complex, receptor, csv -> [ complex, receptor, csv.splitCsv(header:true, strip:true) ] }
                              .map{ complex, receptor, row -> [ complex, receptor, row.Tool, row.Complex, row.Pocket, row.Rank, row.'lDDT-PLI', row.'lDDT-LP', row.BiSyRMSD, row.Reference_Ligand ] }
                              .transpose()
                              .collectFile() { item -> [ "${item[1]}____${item[0]}_${item[4]}_${item[5]}_edmdock_score_summary.csv", item[2] + "," + item[3] + "," + item[4] + "," + item[5] +  "," + item[6] + "," + item[7] + "," + item[8] + "," + item[9] + "\n" ]}
                              .map{ [ it.simpleName.split('____')[0], it.simpleName.split('_')[-5], it.simpleName.split('_pocket')[0].split('____')[1], it ] }
                              .set{ edm_scores_for_coordinates }    // receptor, pocket, complex, edmdock_score_csv

                edm_scores_for_coordinates.combine( predicted_coordinates, by: [0, 1] )        // receptor, pocket, complex, edmdock_score_csv, x, y, z
                                          .set{ edm_scores_for_coordinates_p2rank }

                edm_scores_for_coordinates.combine( defined_coordinates, by:  [0, 1, 2] )
                                          .set{ edm_scores_for_coordinates_defined }

                edm_scores_for_coordinates_p2rank.concat( edm_scores_for_coordinates_defined )
                            .map{ [ it[3].simpleName.split('____')[1], it[0], it[1], it[2], it[3], it[4], it[5], it[6] ] }
                            .map{ file_name, receptor, pocket, complex, csv, x, y, z -> [ file_name, receptor, pocket, complex, csv.splitCsv(strip:true), x, y, z] }
                            .map{ file_name, receptor, pocket, complex, row, x, y, z -> [ file_name, receptor, pocket, complex, row[0][0], row[0][1], row[0][2], row[0][3], row[0][4], row[0][5], row[0][6], row[0][7], x, y, z ] }
                            .transpose()
                            .collectFile() { item -> [ "${item[0]}.csv", item[4] + "," + item[5] + "," + item[6] + "," + item[7] + "," + item[8] +  "," + item[9] + "," + item[10] + "," + item[11] + "," + item[12] + "," + item[13] + "," + item[14] + '\n'] }
                            .collect()
                            .set{ edmdock_scores_for_summary}
             }


            // create ligand score summary file
            overall_ost_scores = ost_score_summary( tb_scores_for_summary.ifEmpty( [] ).combine(
                                                   dd_scores_for_summary.ifEmpty( [] ).combine(
                                                   vina_scores_for_summary.ifEmpty( [] ).combine(
                                                   smina_scores_for_summary.ifEmpty( [] ).combine(
                                                   gnina_scores_for_summary.ifEmpty( [] ).combine(
                                                   edmdock_scores_for_summary.ifEmpty( [] )))))))
        }
        else {
            overall_ost_scores = Channel.empty()
        }


        /*
        * combine tool scores with ost scores summary file
        */

        // edmdock has no own score --> if only edmdock run, no tool score summary needed
        if ( params.tools != /edmdock/ ) {
            overall_ost_scores.ifEmpty( [] ).combine( for_tb_tool_scores.ifEmpty( [] )
                                            .combine( for_dd_tool_scores.ifEmpty( [] )
                                            .combine( for_vina_tool_scores.ifEmpty( [] )
                                            .combine( for_smina_tool_scores.ifEmpty( [] )
                                            .combine( for_gnina_tool_scores.ifEmpty( [] )
                                            .combine( for_edmdock_tool_scores.ifEmpty( [] ) ))))))
                                            .flatten()
                                            .branch{
                                                for_linking: it.name == 'ligand_score_summary.csv' || it.name == 'dd_file_names.txt'
                                                other: true
                                                }
                                            .set{ all_scores_input }

            all_scores = combine_all_scores( all_scores_input.for_linking.collect() )
        }
        else {
            overall_ost_scores.ifEmpty( [] ).map{ [ [it], "True" ] }.map{ it[1] }
                                            .branch{
                                                ready: it == "True"
                                                other: true
                                            }
                                            .set{ all_scores }
        }

        /*
        * error reports
        */

        ignored_tasks = catch_ignored_tasks( all_scores.ready )
        error_and_problems_summary( ignored_tasks.concat( ligand_prep_log, box_size_failed, p2rank_no_pocket, dd_problems ).toList().map{ [ it, "${params.P2RANK}" ] } )
    }
    else {
        binding_pockets.pockets.map{ receptor, csv -> [ receptor, csv.splitCsv(header:true, strip:true) ] }
                               .branch{
                                    p2rank_ok: it[1] != []
                                    other: true
                               }
                               .set{ p2rank_success }

        p2rank_success.other.collectFile( storeDir: 'preprocessing/p2rank' ) {item -> [ "p2rank_no_pockets_found.csv", item[0] + "\n"] }
                            .set{ p2rank_no_pocket }

        ignored_tasks = catch_ignored_tasks( binding_pockets.pockets.collect().map{ [ [it], "True" ] }.map{ it[1] } )
        error_and_problems_summary( ignored_tasks.concat( p2rank_no_pocket ).toList().map{ [ it, "${params.P2RANK}" ] } )
    }

}
