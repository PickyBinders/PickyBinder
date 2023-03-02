/*
*  add_Hs_to_receptor module 
*/

params.OUTPUT = "$launchDir/data/receptor_Hs"

process add_Hs_to_receptor {
    publishDir "$params.OUTPUT", mode: 'copy'
    container "/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/ADFRsuite.img"
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_file)

    output:
    tuple val (receptor), path ("*_Hs.pdb"), emit: pdb_Hs

    script:
    """
    reduce ${pdb_file} > ${receptor}_Hs.pdb
    """
}
