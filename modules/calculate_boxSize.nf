/*
*  calculate_boxSize module 
*/

params.OUTPUT = "$launchDir/boxes"

process calculate_boxSize {
    publishDir("$launchDir/boxes/box_size", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina_meeko'
    tag { ligand }

    input:
    tuple val (ligand), path (ligand_sdf)

    output:
    tuple path ("${ligand}_boxSize.txt"), val (ligand), emit: box_size

    script:
    """
    calculate_boxSize.py ${ligand_sdf} ${ligand} ${params.autobox_add} > ${ligand}_boxSize.txt
    """
}
