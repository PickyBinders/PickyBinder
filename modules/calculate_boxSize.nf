/*
*  calculate_boxSize module 
*/

params.OUTPUT = "$launchDir/boxes/box_size"

process calculate_boxSize {
    publishDir "$params.OUTPUT", pattern: "*_boxSize.txt", mode: 'copy'
    publishDir "$params.OUTPUT/log_files", pattern: "*_log.txt", mode: 'copy'
    conda "${params.meeko_conda}"
    tag { ligand }

    input:
    tuple val (ligand), path (ligand_sdf)

    output:
    tuple path ("${ligand}_boxSize.txt"), val (ligand), optional: true, emit: size
    path ("${ligand}_log.txt"), optional: true

    script:
    """
    calculate_boxSize.py ${ligand_sdf} ${ligand} ${params.autobox_add} &> ${ligand}_log.txt
    """
}
