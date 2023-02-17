/*
*  p2rank module 
*/

//params.OUTPUT = "$launchDir/p2rank"

process p2rank {
    publishDir("$launchDir/p2rank/${receptor_chain}", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'
    tag { receptor_chain }

    input:
    tuple val (receptor_chain), path (pdb_Hs)

    output:
    tuple val (receptor_chain), path ("*_predictions.csv"), emit: pockets
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
