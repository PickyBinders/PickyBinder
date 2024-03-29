/*
*  p2rank module 
*/

params.OUTPUT = "$launchDir/preprocessing/p2rank"

process p2rank {
    publishDir "$params.OUTPUT/${receptor}", mode: 'copy'
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_Hs), val (alphafold)

    output:
    tuple val (receptor), path ("*_predictions.csv"), emit: pockets
    path ("*_residues.csv")
    path ("visualizations/")
    path ("params.txt")
    path ("run.log")

    script:
    """
    ${params.load_java}

    if test "${alphafold}" = "yes"
    then
        ${params.p2rank_tool}/prank predict -f ${pdb_Hs} -o ./ -c alphafold
    else
        ${params.p2rank_tool}/prank predict -f ${pdb_Hs} -o ./
    fi
    """
}
