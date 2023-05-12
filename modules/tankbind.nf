/*
*  tankbind module
*/

params.OUTPUT = "$launchDir/tankbind"

process tankbind {
    publishDir "$params.OUTPUT/${complex}", pattern: "tankbind_predictions/*", mode: 'copy'
    publishDir "$params.OUTPUT/${complex}", pattern: "*_tankbind.csv", mode: 'copy'
    publishDir "$params.OUTPUT/${complex}", mode: 'copy', saveAs: { filename -> if (filename == ".command.log") "tankbind.log"}
    container "${params.tankbind_sing}"
    containerOptions "-B ${params.tankbind_scripts}"
    tag { complex }
    
    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_Hs), path (p2rank_prediction), path (ligand_sdf)
    
    output:
    tuple val (complex), path ("${complex}_tankbind.csv"), emit: affinities
    tuple val (complex), path ("tankbind_predictions/*"), emit: sdfs
    path (".command.log"), emit: tankbind_log
    
    script:
    """
    tankbind_prediction.py ${params.tankbind_scripts} ${complex} ${pdb_Hs} ${ligand_sdf} ${p2rank_prediction} ${params.naming}
    """
}
