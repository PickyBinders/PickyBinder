/*
*  gnina module
*/

params.OUTPUT = "$launchDir/gnina"


process gnina {
    publishDir("$launchDir/gnina/gnina_predictions/${complex}/${pocket_nr}", mode: 'copy')
    container '/scicore/home/schwede/leeman0000/singularity/nmaus-gnina.img'
    tag { complex }
    
    input:
    tuple val (complex), val (ligand), val (receptor_chain), val (pocket_nr), path (receptor_pdbqt), path (ligand_pdbqt), path (vina_box)
    
    output:
    tuple val (complex), val (receptor_chain), val (pocket_nr), path ("${complex}_${pocket_nr}_gnina.pdbqt"), emit: gnina_result
    path ("${complex}_${pocket_nr}_gnina.log"), emit: gnina_log
    
    script:
    """
    gnina -r ${receptor_pdbqt} -l ${ligand_pdbqt} --config ${vina_box} \
          -o ${complex}_${pocket_nr}_gnina.pdbqt \
          --log ${complex}_${pocket_nr}_gnina.log \
          --exhaustiveness=32 --seed 160490
    """
}