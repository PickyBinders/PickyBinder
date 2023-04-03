/*
*  smina module
*/

params.OUTPUT = "$launchDir/smina"

process smina {
    publishDir "$params.OUTPUT/smina/${complex}/${pocket_nr}", mode: 'copy'
    conda "/scicore/home/schwede/leeman0000/miniconda3/envs/smina"
    tag { complex }
    
    input:
    tuple val (complex), val (ligand), val (receptor), val (pocket_nr), path (receptor_pdbqt), path (ligand_pdbqt), path (vina_box)
    
    output:
    tuple val (complex), val (receptor), val (pocket_nr), path ("${complex}_${pocket_nr}_smina.pdbqt"), emit: smina_result
    path ("${complex}_${pocket_nr}_smina.log"), emit: smina_log
    
    script:
    """
    smina -r ${receptor_pdbqt} -l ${ligand_pdbqt} --config ${vina_box} \
          -o ${complex}_${pocket_nr}_smina.pdbqt \
          --log ${complex}_${pocket_nr}_smina.log \
          --exhaustiveness=64 --seed 160490 --cpu ${task.cpus}
    """
}


process smina_sdf {
    publishDir "$params.OUTPUT/smina/${complex}/${pocket_nr}", mode: 'copy'
    conda "/scicore/home/schwede/leeman0000/miniconda3/envs/smina"
    tag { complex }

    input:
    tuple val (complex), val (ligand), val (receptor), val (pocket_nr), path (receptor_pdb), path (ligand_sdf), path (vina_box)

    output:
    tuple val (complex), val (receptor), val (pocket_nr), path ("${complex}_${pocket_nr}_smina_*.sdf"), emit: smina_sdf
    path ("${complex}_${pocket_nr}_smina.log"), emit: smina_log

    script:
    """
    smina -r ${receptor_pdb} -l ${ligand_sdf} --config ${vina_box} \
          -o ${complex}_${pocket_nr}_smina.sdf \
          --log ${complex}_${pocket_nr}_smina.log \
          --exhaustiveness=64 --seed 160490 --cpu ${task.cpus}

    split_pat='\$\$\$\$'
    csplit --elide-empty-files --prefix="${complex}_${pocket_nr}_smina_" --suffix-format="%d.sdf"  <(echo \$split_pat; cat "${complex}_${pocket_nr}_smina.sdf") '/\$\$\$\$/+1' "{*}"
    for line in ${complex}_${pocket_nr}_smina_*; do if test \$(wc -l < \$line) -eq 1; then rm \$line; fi; done
    """
}
