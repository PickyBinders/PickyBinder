/*
*  docking_box module 
*/

params.OUTPUT = "$launchDir/boxes"


process docking_boxes_predicted_pockets {
    publishDir "$params.OUTPUT/${complex}", mode: 'copy'
    container "${params.tankbind_sing}"
    containerOptions "-B ${params.tankbind_scripts}"
    tag { complex }

    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_Hs), path (p2rank_predictions), path (box_size)

    output:
    path ("*.box"), emit: box_per_pocket

    script:
    """
    define_docking_boxes.py ${params.tankbind_scripts} ${complex} ${pdb_Hs} ${box_size} ${p2rank_predictions} false
    """
}


process docking_box_defined_BS {
    publishDir "$params.OUTPUT/${complex}", mode: 'copy'
    container "${params.tankbind_sing}"
    containerOptions "-B ${params.tankbind_scripts}"
    tag { complex }

    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_Hs), val (coordinates), path (box_size)

    output:
    path ("*.box"), emit: boxes_bs_wp

    script:
    """
    define_docking_boxes.py ${params.tankbind_scripts} ${complex} ${pdb_Hs} ${box_size} false ${coordinates}
    """
}