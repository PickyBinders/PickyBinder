/*
*  p2rank module 
*/

params.OUTPUT = "$launchDir/p2rank"

process p2rank {
    publishDir "$params.OUTPUT/${receptor}", mode: 'copy'
    conda "/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd"
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_Hs)

    output:
    tuple val (receptor), path ("*_predictions.csv"), emit: pockets
    path ("*_residues.csv")
    path ("visualizations/")
    path ("params.txt")
    path ("run.log")

    script:
    """
    ml load Java/13.0.2
    
    ${params.p2rank_tool}/prank predict -f ${pdb_Hs} -o ./
    """
}
