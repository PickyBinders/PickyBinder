/*
*  tankbind module
*/

params.OUTPUT = "$launchDir/predictions/tankbind"

process tankbind {
    publishDir "$params.OUTPUT/${complex}", pattern: "tankbind_predictions/*", mode: 'copy'
    publishDir "$params.OUTPUT/${complex}", pattern: "*_tankbind.csv", mode: 'copy'
    publishDir "$params.OUTPUT/${complex}", mode: 'copy', pattern: "tankbind_*.log"
    container "${params.tankbind_sing}"
    containerOptions "-B ${params.tankbind_scripts}"
    tag { complex }
    
    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_Hs), path (p2rank_prediction), path (ligand_sdf)
    
    output:
    tuple val (complex), val (receptor), path ("${complex}_tankbind.csv"), emit: affinities
    tuple val (complex), path ("tankbind_predictions/*"), emit: sdfs
    path ("tankbind_*.log"), emit: tankbind_log
    
    script:
    """
    tankbind_prediction.py ${params.tankbind_scripts} ${complex} ${receptor} ${pdb_Hs} ${ligand} ${ligand_sdf} ${p2rank_prediction}

    ln -s .command.log tankbind_${complex}.log
    """
}
