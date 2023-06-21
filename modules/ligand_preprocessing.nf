/*
*  ligand_preprocessing module
*/

params.OUTPUT = "$launchDir/preprocessing/ligands"

process ligand_preprocessing {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.{sdf,txt}"
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "ligand_preprocessing_*.log"
    conda "${params.meeko_conda}"

    input:
    path (ref_sdf_files)
    path(mol_files)

    output:
    path ("*.sdf"), emit: sdf_files
    path ("problems_ligand_prep_*.txt"), emit: lig_prep_problems
    path ("ligand_preprocessing_*.log"), emit: ligand_prep_log

    script:
    """
    ligand_preprocessing.py ${params.naming} ${ref_sdf_files}

    time_date=\$(date +"%y-%m-%d-%T")
    ln -s .command.log ligand_preprocessing_\${time_date}.log
    mv problems_ligand_prep.txt problems_ligand_prep_\${time_date}.txt
    """
}