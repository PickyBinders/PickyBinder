/*
*  vina module
*/

params.OUTPUT = "$launchDir/vina"

process vina_prepare_receptor {
    publishDir(params.OUTPUT, mode: 'copy')
    container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/ADFRsuite.img'
    tag { sample_name }

    input:
    tuple val (sample_name), path (receptor_pdb)
    
    output:
    tuple val (sample_name), path ("${sample_name}.pdbqt"), emit: preped_receptor
    
    script:
    """
    prepare_receptor -r ${receptor_pdb} -o ${sample_name}.pdbqt
    """
}


process vina_prepare_ligand {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { sample_name }
    
    input:
    tuple val (sample_name), path (ligand_sdf)
    
    output:
    tuple val (sample_name), path ("${sample_name}.pdbqt"), emit: preped_ligand
    
    script:
    """
    mk_prepare_ligand.py -i ${ligand_sdf} -o ${sample_name}.pdbqt
    """
}


process vina_box {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/mgltools'
    tag { sample_name }
    
    input:
    path (receptor_pdbqt)
    tuple val (sample_name), path (ligand_pdbqt)
    
    output:
    path ("${sample_name}.gpf"), emit: grid_file
    path ("${sample_name}_box.txt"), emit: box
    
    script:
    """
    $baseDir/bin/vina_scripts/prepare_gpf_nextflow.py -l ${ligand_pdbqt} -r ${receptor_pdbqt} -y
    
    echo -e "center_x = \$(grep gridcenter ${sample_name}.gpf | cut -d " " -f2) \
            \ncenter_y = \$(grep gridcenter ${sample_name}.gpf | cut -d " " -f3) \
            \ncenter_z = \$(grep gridcenter ${sample_name}.gpf | cut -d " " -f4) \
            \nsize_x = 20.0 \nsize_y = 20.0 \nsize_z = 20.0 \n" >> ${sample_name}_box.txt
    """
}


process vina {
    publishDir(params.OUTPUT, mode: 'copy')
    //container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/AutoDockVina.img' 
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina'
    tag { sample_name }
    
    input:
    path (receptor_pdbqt)
    path (ligand_pdbqt)
    path (vina_box)
    
    output:
    path ("*vina_out.pdbqt"), emit: vina_result
    path ("vina.log"), emit: vina_log
    
    script:
    """
    /scicore/home/schwede/leeman0000/tools/vina/vina_1.2.3_linux_x86_64 \
         --receptor ${receptor_pdbqt} --ligand ${ligand_pdbqt} \
         --config ${vina_box} \
         --exhaustiveness=32 --out vina_out.pdbqt \
         > vina.log
    """
}