/*
*  docking_box module 
*/

params.OUTPUT = "$launchDir/boxes"

process docking_box {
    publishDir("$launchDir/boxes/${complex}", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina_meeko'
    tag { ligand }

    input:
    tuple val (ligand), val (receptor_chain), val (complex), path (p2rank_predictions), path (box_size)

    output:
    path ("*.box"), emit: box_per_pocket

    script:
    """
    define_box_per_pocket.py ${complex} ${p2rank_predictions} ${box_size}
    """
}
