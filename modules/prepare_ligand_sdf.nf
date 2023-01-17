/*
*  prepare_ligand_sdf module 
*/

params.OUTPUT = "$launchDir/data/ligands"

process prepare_ligand_sdf {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'

    input:
    path (ref_sdf_files)

    output:
    path ("*.sdf"), emit: sdf_files

    script:
    """
    prepare_ligand_sdf.py ${ref_sdf_files}
    """
}
