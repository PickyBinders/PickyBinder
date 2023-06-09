/*
*  prepare_ligand_sdf module 
*/

params.OUTPUT = "$launchDir/preprocessing/ligands"

process prepare_ligand_sdf {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.sdf"
    publishDir "$params.OUTPUT", mode: 'copy', saveAs: { filename -> if (filename == ".command.log") "ligand_preparation.log"}
    conda "${params.meeko_conda}"

    input:
    path (ref_sdf_files)
    path(mol_files)

    output:
    path ("*.sdf"), emit: sdf_files
    path (".command.log"), emit: ligand_prep_log

    script:
    """
    prepare_ligand_sdf.py ${ref_sdf_files}
    """
}


process ligand_preprocessing {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.{sdf,txt}"
    //publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.txt"
    publishDir "$params.OUTPUT", mode: 'copy', saveAs: { filename -> if (filename == ".command.log") "ligand_preprocessing.log"}
    conda "${params.meeko_conda}"

    input:
    path (ref_sdf_files)
    path(mol_files)

    output:
    path ("*.sdf"), emit: sdf_files
    path ("problems_ligand_prep.txt"), emit: lig_prep_problems
    path (".command.log"), emit: ligand_prep_log

    script:
    """
    ligand_preprocessing.py ${ref_sdf_files}
    """
}