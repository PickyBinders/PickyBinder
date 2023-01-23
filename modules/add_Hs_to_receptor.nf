/*
*  add_Hs_to_receptor module 
*/

params.OUTPUT = "$launchDir/data/hydrated_receptors"

process add_Hs_to_receptor {
    publishDir(params.OUTPUT, mode: 'copy')
    container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/ADFRsuite.img'
    tag { receptor_chain }

    input:
    tuple val (receptor_chain), path (pdb_file)

    output:
    tuple val (receptor_chain), path ("*_Hs.pdb"), emit: pdb_H_files

    script:
    """
    reduce ${pdb_file} > ${receptor_chain}_Hs.pdb
    """
}
