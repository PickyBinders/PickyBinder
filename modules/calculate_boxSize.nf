/*
*  calculate_boxSize module 
*/

params.OUTPUT = "$launchDir/boxes"

process calculate_boxSize {
    publishDir("$launchDir/boxes/box_size", pattern: "*_boxSize.txt", mode: 'copy')
    publishDir("$launchDir/boxes/box_size/failed", pattern: "*_log.txt", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina_meeko'
    tag { ligand }

    input:
    tuple val (ligand), path (ligand_sdf)

    output:
    tuple path ("${ligand}_boxSize.txt"), val (ligand), optional: true, emit: size
    path ("${ligand}_log.txt"), optional: true

    script:
    """
    calculate_boxSize.py ${ligand_sdf} ${ligand} ${params.autobox_add} &> ${ligand}_log.txt

    if [ -s ${ligand}_boxSize.txt ]; then rm ${ligand}_log.txt; fi
    """
}
