/*
*  vina module
*/

params.OUTPUT = "$launchDir/vina"


process vina_prepare_ligand {
    //publishDir(params.OUTPUT, mode: 'copy', pattern: "*.pdbqt")
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), val (recep_chain), path (sdf_file)
    path (pdb_files)
    
    output:
    tuple val (rec_lig), path ("${rec_lig}.pdbqt"), val (recep_chain), path ("${recep_chain}.pdb"), emit: preped_ligand
    
    script:
    """
    mk_prepare_ligand.py -i ${sdf_file} -o ${rec_lig}.pdbqt
    
    cp ${pdb_files}/${recep_chain}.pdb ${recep_chain}.pdb
    """
}


process vina_prepare_receptor {
    //publishDir(params.OUTPUT, mode: 'copy')
    container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/ADFRsuite.img'
    tag { rec_lig }

    input:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdb)
    
    output:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path ("${recep_chain}.pdbqt"), emit: preped_receptor
    
    script:
    """
    prepare_receptor -r ${receptor_pdb} -o ${recep_chain}.pdbqt
    """
}


process vina_box {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/mgltools'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt)
    
    output:
    path ("${rec_lig}.gpf"), emit: grid_file
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt), path ("${rec_lig}_box.txt"), emit: pdbqtFiles_box
    
    script:
    """
    $baseDir/bin/vina_scripts/prepare_gpf_nextflow.py -l ${ligand_pdbqt} -r ${receptor_pdbqt} -y -o ${rec_lig}.gpf
    
    echo -e "center_x = \$(grep gridcenter ${rec_lig}.gpf | cut -d " " -f2) \
            \ncenter_y = \$(grep gridcenter ${rec_lig}.gpf | cut -d " " -f3) \
            \ncenter_z = \$(grep gridcenter ${rec_lig}.gpf | cut -d " " -f4) \
            \nsize_x = 20.0 \nsize_y = 20.0 \nsize_z = 20.0 \n" >> ${rec_lig}_box.txt
    """
}


process vina {
    publishDir(params.OUTPUT, mode: 'copy')
    //container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/AutoDockVina.img' 
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt), path (vina_box)
    
    output:
    tuple val (rec_lig), val (recep_chain), path ("${rec_lig}_vina.pdbqt"), emit: vina_result
    path ("${rec_lig}_vina.log"), emit: vina_log
    
    script:
    """
    /scicore/home/schwede/leeman0000/tools/vina/vina_1.2.3_linux_x86_64 \
         --receptor ${receptor_pdbqt} --ligand ${ligand_pdbqt} \
         --config ${vina_box} \
         --exhaustiveness=32 --out ${rec_lig}_vina.pdbqt \
         > ${rec_lig}_vina.log
    """
}