/*
*  docking_box module 
*/

params.OUTPUT = "$launchDir/boxes"

process docking_box {
    publishDir "$params.OUTPUT/${complex}", mode: 'copy'
    container "/scicore/home/schwede/leeman0000/singularity/qizhipei-tankbind_py38.img"
    containerOptions "-B ${params.tankbind_scripts}"
    tag { ligand }

    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_Hs), path (p2rank_predictions), path (box_size)

    output:
    path ("*.box"), emit: box_per_pocket

    script:
    """
    define_box_per_pocket.py ${params.tankbind_scripts} ${complex} ${p2rank_predictions} ${box_size} ${pdb_Hs}
    """
}
