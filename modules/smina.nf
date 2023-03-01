/*
*  smina module
*/

params.OUTPUT = "$launchDir/smina"


process smina {
    publishDir("$launchDir/smina/smina_predictions/${complex}/${pocket_nr}", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/smina'
    tag { complex }
    
    input:
    tuple val (complex), val (ligand), val (receptor_chain), val (pocket_nr), path (receptor_pdbqt), path (ligand_pdbqt), path (vina_box)
    
    output:
    tuple val (complex), val (receptor_chain), val (pocket_nr), path ("${complex}_${pocket_nr}_smina.pdbqt"), emit: smina_result
    path ("${complex}_${pocket_nr}_smina.log"), emit: smina_log
    
    script:
    """
    smina -r ${receptor_pdbqt} -l ${ligand_pdbqt} --config ${vina_box} \
          -o ${complex}_${pocket_nr}_smina.pdbqt \
          --log ${complex}_${pocket_nr}_smina.log \
          --exhaustiveness=32 --seed 160490 --cpu ${task.cpus}
    """
}