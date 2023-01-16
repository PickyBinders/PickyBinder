/*
*  vina module
*/

params.OUTPUT = "$launchDir/vina"


process vina_prepare_ligand {
    //publishDir(params.OUTPUT, mode: 'copy', pattern: "*.pdbqt")
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { rec_lig }
    
    input:
    tuple val (lig), val (rec_lig), val (recep_chain) 
    path (sdf_file)
    path (pdb_files)
    
    output:
    tuple val (rec_lig), path ("${rec_lig}.pdbqt"), val (recep_chain), path ("${recep_chain}.pdb"), emit: preped_ligand
    
    script:
    """
    mk_prepare_ligand.py -i ${lig}.sdf -o ${rec_lig}.pdbqt
    
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
    reduce ${receptor_pdb} > ${recep_chain}_H.pdb
    prepare_receptor -r ${recep_chain}_H.pdb -o ${recep_chain}.pdbqt
    """
}


process vina_box {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt)
    path (receptors)
    
    output:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt), path ("${rec_lig}_box.txt"), emit: pdbqtFiles_box
    
    script:
    """
    calculate_box_for_vina.py ${recep_chain} ${rec_lig}
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
    ${params.vina_tool} --receptor ${receptor_pdbqt} --ligand ${ligand_pdbqt} \
                 --config ${vina_box} \
                 --exhaustiveness=32 --out ${rec_lig}_vina.pdbqt \
                 > ${rec_lig}_vina.log
    """
}    


process vina_pdbtqToSdf {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), val (recep_chain), path (vina_pdbqt)
    
    output:
    tuple val (rec_lig), val (recep_chain), path ("${rec_lig}_vina*.sdf")
    
    shell:
    '''
    mk_export.py !{vina_pdbqt} -o !{rec_lig}_vina.sdf
    
    csplit --prefix="!{rec_lig}_vina_" --suffix-format="%d.sdf"  <(echo "\\$\\$\\$\\$"; cat !{rec_lig}_vina.sdf) '/$$$$/+1' "{*}"
    for line in !{rec_lig}_vina_*; do if test $(wc -l < $line) -eq 1; then rm $line; fi; done
    '''
}


/* #!/usr/bin/env python
*    
*    from meeko import PDBQTMolecule
*    from rdkit import Chem
*    
*    pdbqt_mol = PDBQTMolecule.from_file('${vina_pdbqt}', is_dlg=False, skip_typing=True)
*    for i, pose in enumerate(pdbqt_mol):
*        rdkit_mol = pose.export_rdkit_mol()
*        affinity = pose.score
*        with Chem.SDWriter(f"${rec_lig}_vina_{i}.sdf") as w:
*            w.write(rdkit_mol)
*/


process vina_all {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina_meeko'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), path (ligand_pdbqt), val (recep_chain), path (receptor_pdbqt)
    path (receptors)
    
    output:
    tuple val (rec_lig), val (recep_chain), path ("${rec_lig}_*.sdf"), emit: vina_result
    
    script:
    """
    run_vina.py ${recep_chain} ${rec_lig} ${receptor_pdbqt} ${ligand_pdbqt}
    """

}


