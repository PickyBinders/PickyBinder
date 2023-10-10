/*
*  ligand_preprocessing module
*/

params.OUTPUT = "$launchDir/preprocessing/ligands"

process ligand_preprocessing_single {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "*.{sdf}"
    publishDir "$params.OUTPUT/log_files", mode: 'copy', pattern: "*_ligand_preprocessing.log"
    conda "${params.meeko_conda}"
    tag { ligand }

    input:
    tuple val(ligand), path (ref_sdf_file), path(mol_file)

    output:
    path ("*.sdf"), optional: true, emit: sdf_files
    path ("*_ligand_preprocessing.log"), emit: ligand_prep_log

    script:
    """
    ligand_preprocessing.py ${ref_sdf_file} ${ligand}

    ln -s .command.log ${ligand}_ligand_preprocessing.log
    """
}


process ligand_preprocessing_log {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "ligand_preprocessing_*.log"

    input:
    path (all_log_files)

    output:
    path ("ligand_preprocessing_*.log")

    script:
    """
    time_date=\$(date +"%y-%m-%d-%T")

    echo -e "**************\nLigand preparation failed:" > log.txt
    grep -h -A1 'Ligand preparation failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Ligand preparation failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nUncharging failed:" >> log.txt
    grep -h -A1 'Uncharging failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Uncharging failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nProtonation failed:" >> log.txt
    grep -h -A1 'Protonation failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Protonation failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nEmbedding failed:" >> log.txt
    grep -h -A1 'Embedding failed:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Embedding failed:' | tr -s '\n' >> log.txt
    echo -e "\n**************\nWarnings:" >> log.txt
    grep -h -A1 'Warnings:' *_ligand_preprocessing.log --no-group-separator | grep -v 'Warnings:' | tr -s '\n' >> log.txt
    echo -e "\n" >> log.txt
    mv log.txt ligand_preprocessing_\${time_date}.log
    """
}
