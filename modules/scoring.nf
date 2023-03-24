/*
*  scoring module 
*/

params.OUTPUT = "$launchDir/scores"

process ost_scoring_models {
    publishDir "$params.OUTPUT/${complex}/${tool_name}", mode: 'copy'
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor), path (model_receptor), path (ref_ligand), path (modeled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score
    tuple val (complex), path ("*.csv"), emit: summary

    script:
    """
    source /scicore/home/schwede/leeman0000/activate-ost-develop

    for model in ${modeled_ligands}
        do
        ost compare-ligand-structures -m ${model_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} -o \${model%.sdf}.json --lddt-pli --rmsd
        done

    combine_ost_scores.py ${tool_name}
    """
}


process ost_scoring {
    publishDir "$params.OUTPUT/${complex}/${tool_name}", mode: 'copy'
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor), path (ref_ligand), path (modeled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score
    tuple val (complex), path ("*.csv"), emit: summary

    script:
    """
    source /scicore/home/schwede/leeman0000/activate-ost-develop

    for model in ${modeled_ligands}
        do
        ost compare-ligand-structures -m ${ref_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} -o \${model%.sdf}.json --lddt-pli --rmsd
        done

    combine_ost_scores.py ${tool_name}
    """
}
