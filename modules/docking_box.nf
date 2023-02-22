/*
*  docking_box module 
*/

params.OUTPUT = "$launchDir/boxes"

process docking_box {
    publishDir("$launchDir/boxes/${complex}", mode: 'copy')
    //conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina_meeko'
    container '/scicore/home/schwede/leeman0000/singularity/qizhipei-tankbind_py38.img'
    tag { ligand }
    containerOptions "-B ${params.tankbind_scripts}"

    input:
    tuple val (ligand), val (receptor_chain), val (complex), path (pdb_Hs), path (p2rank_predictions), path (box_size)

    output:
    path ("*.box"), emit: box_per_pocket

    script:
    """
    define_box_per_pocket.py ${params.tankbind_scripts} ${complex} ${p2rank_predictions} ${box_size} ${pdb_Hs}
    """
}
