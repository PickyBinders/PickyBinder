/*
*  scoring module 
*/

params.OUTPUT = "$launchDir/scores"

process ost_scoring {
    publishDir "$params.OUTPUT/${complex}/${tool_name}", mode: 'copy'
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (receptor), path (ref_ligand), path (modeled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score

    script:
    """
    source /scicore/home/schwede/leeman0000/activate-ost-develop

    for model in ${modeled_ligands}
        do
        ost compare-ligand-structures -m ${receptor} -ml \${model} -r ${receptor} -rl ${ref_ligand} -o \${model}.json --lddt-pli --rmsd
        done
    """
}
