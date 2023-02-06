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
    publishDir("$launchDir/vina/vina_predictions/${rec_lig}", mode: 'copy')
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
                 > ${rec_lig}_vina.log --seed 160490
    """
}    


process vina_pdbtqToSdf {
    publishDir("$launchDir/vina/vina_predictions/${rec_lig}", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { rec_lig }
    
    input:
    tuple val (rec_lig), val (recep_chain), path (vina_pdbqt)
    
    output:
    tuple val (rec_lig), val (recep_chain), path ("${rec_lig}_vina*.sdf")
    
    script:
    """
    csplit --elide-empty-files --prefix="${rec_lig}_vina_" --suffix-format="%d.pdbqt"  <(echo "ENDMDL"; cat ${vina_pdbqt}) '/ENDMDL/+1' "{*}"
    for line in ${rec_lig}_vina_*; do if test \$(wc -l < \$line) -eq 1; then rm \$line; fi; done
     
    for file in ${rec_lig}_vina_*.pdbqt; do mk_export.py \${file} -o \${file%.pdbqt}.sdf; done
    """
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


process vina_prepare_ligand2 {
    publishDir("$launchDir/vina/prepared_ligands", mode: 'copy', pattern: "*.pdbqt")
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    
    input:
    path (sdf_files)
    
    output:
    path ("*.pdbqt"), emit: preped_ligand
    
    script:
    """
    for file in *.sdf; do mk_prepare_ligand.py -i \${file} -o \${file%.sdf}.pdbqt; done
    """
}


process vina_prepare_receptor2 {
    publishDir("$launchDir/vina/prepared_receptors", mode: 'copy')
    container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/ADFRsuite.img'
    tag { receptor_chain }

    input:
    tuple val(receptor_chain), path (pdb_H_files)
    
    output:
    tuple val(receptor_chain), path ("*.pdbqt"), emit: preped_receptor
    
    script:
    """
    prepare_receptor -r ${pdb_H_files} -o ${receptor_chain}.pdbqt
    """
}


process vina_box2 {
    publishDir("$launchDir/vina/box_files", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'
    tag { receptor_chain }
    
    input:
    tuple val (receptor_chain), path (pdb_Hs)

    output:
    tuple val (receptor_chain), path ("*_box.txt"), emit: box
    
    script:
    """
    calculate_box_for_vina2.py ${receptor_chain} ${pdb_Hs}
    """
}


process vina2 {
    publishDir("$launchDir/vina/vina_predictions/${complex}", mode: 'copy')
    //container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/AutoDockVina.img' 
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina'
    tag { complex }
    
    input:
    tuple val (ligand), val (receptor_chain), val (complex), path (receptor_pdbqt), path (vina_box), path (ligand_pdbqt)
    
    output:
    tuple val (complex), val (receptor_chain), path ("${complex}_vina.pdbqt"), emit: vina_result
    path ("${complex}_vina.log"), emit: vina_log
    
    script:
    """
    ${params.vina_tool} --receptor ${receptor_pdbqt} --ligand ${ligand_pdbqt} \
                 --config ${vina_box} \
                 --exhaustiveness=32 --out ${complex}_vina.pdbqt \
                 > ${complex}_vina.log --seed 160490
    """
}


process vina3 {
    publishDir("$launchDir/vina/vina_predictions/${complex}/${pocket_nr}", mode: 'copy')
    //container '/scicore/home/schwede/zohixe92/CAMEO/CAMEO_predictors/BaselineCM_AutoDockVina/container_scripts/AutoDockVina.img' 
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/vina'
    tag { complex }
    
    input:
    tuple val (complex), val (ligand), val (receptor_chain), val (pocket_nr), path (receptor_pdbqt), path (ligand_pdbqt), path (vina_box)
    
    output:
    tuple val (complex), val (receptor_chain), val (pocket_nr), path ("${complex}_${pocket_nr}_vina.pdbqt"), emit: vina_result
    path ("${complex}_${pocket_nr}_vina.log"), emit: vina_log
    
    script:
    """
    ${params.vina_tool} --receptor ${receptor_pdbqt} --ligand ${ligand_pdbqt} \
                 --config ${vina_box} --out ${complex}_${pocket_nr}_vina.pdbqt \
                 --exhaustiveness=32 --seed 160490 \
                 > ${complex}_${pocket_nr}_vina.log
    """
}


process vina_pdbtqToSdf3 {
    publishDir("$launchDir/${tool_name}/${tool_name}_predictions/${complex}/${pocket_nr}", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/meeko'
    tag { complex }
    
    input:
    tuple val (complex), val (receptor_chain), val (pocket_nr), path (pdbqt)
    val (tool_name)
    
    output:
    tuple val (complex), val (receptor_chain), val (pocket_nr), path ("${complex}_${pocket_nr}_${tool_name}*.sdf")
    
    script:
    """
    csplit --elide-empty-files --prefix="${complex}_${pocket_nr}_${tool_name}_" --suffix-format="%d.pdbqt"  <(echo "ENDMDL"; cat ${pdbqt}) '/ENDMDL/+1' "{*}"
    for line in ${complex}_${pocket_nr}_${tool_name}_*; do if test \$(wc -l < \$line) -eq 1; then rm \$line; fi; done
     
    for file in ${complex}_${pocket_nr}_${tool_name}_*.pdbqt; do mk_export.py \${file} -o \${file%.pdbqt}.sdf; done
    """
}


