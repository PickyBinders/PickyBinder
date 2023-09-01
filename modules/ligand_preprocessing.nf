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
    ligand_preprocessing.py ${params.naming} 

    time_date=\$(date +"%y-%m-%d-%T")
    ln -s .command.log ligand_preprocessing_\${time_date}.log
    mv problems_ligand_prep.txt problems_ligand_prep_\${time_date}.txt
    """
}


process ligand_preprocessing_single {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.{sdf}"
    conda "${params.meeko_conda}"
    tag { ligand }

    input:
    tuple val(ligand), path (ref_sdf_file), path(mol_file)

    output:
    path ("*.sdf"), optional: true, emit: sdf_files
    path ("*_ligand_preprocessing.log"), emit: ligand_prep_log

    script:
    """
    ligand_preprocessing.py ${params.naming}

    ln -s .command.log ${ligand}_ligand_preprocessing.log
    """
}

//    path ("problems_ligand_prep_*.txt"), emit: lig_prep_problems
//    mv problems_ligand_prep.txt ${ligand}_problems_ligand_prep.txt
//    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*_ligand_preprocessing.log"

process ligand_preprocessing_log {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "ligand_preprocessing.log"

    input:
    path (all_log_files)

    output:
    path ("ligand_preprocessing.log")

    script:
    """
    echo -e "**************\nUncharging failed:" > log.txt
    grep -h -A1 'Uncharging failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Uncharging failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nProtonation failed:" >> log.txt
    grep -h -A1 'Protonation failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Protonation failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nEmbedding failed:" >> log.txt
    grep -h -A1 'Embedding failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Embedding failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nWarnings:" >> log.txt
    grep -h -A1 'Warnings:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Warnings:' | tr -s '\n' >> log.txt
    echo -e "\n"
    mv log.txt ligand_preprocessing.log
    """
}
