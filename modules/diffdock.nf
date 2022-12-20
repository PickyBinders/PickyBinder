/*
*  diffdock module
*/

params.OUTPUT = "$launchDir"

process diffdock {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/diffdock'
    label 'diffdock'

    input:
    path (protein_ligand_csv)

    output:
    path ("diffdock_predictions/"), emit: diffdock_predictions

    script:
    """
    #!/bin/bash
    
    shopt -s extglob
    ln -s ${params.diffdock_location}/* .
    ln -s ${params.pdb_sdf_files} .
    shopt -u extglob
    
    mkdir data_local
    
    python datasets/esm_embedding_preparation.py --protein_ligand_csv ${protein_ligand_csv} --out_file data_local/prepared_for_esm.fasta
    HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data_local/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok
    
    python -m inference --protein_ligand_csv ${protein_ligand_csv} --out_dir diffdock_predictions \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise
    """
}

process create_csv {
    //publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'
    tag { sample_name }
    
    input:
    tuple val (sample_name), path (sdf_file)
    
    output:
    tuple val (sample_name), path ("${sample_name}_protein_ligand.csv"), emit: protein_ligand_csv
    
    script:
    """
    pdb_name=\$(echo ${sample_name} | cut -d'_' -f1,2)
    
    echo "protein_path,ligand" > ${sample_name}_protein_ligand.csv; \
    echo "${params.pdb_sdf_files}/\${pdb_name}.pdb,${params.pdb_sdf_files}/${sdf_file}" >> ${sample_name}_protein_ligand.csv
    """
} 


process diffdock_single {
    //publishDir(params.OUTPUT, mode: 'copy')
    publishDir("diffdock_predictions", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/diffdock'
    label 'diffdock'
    tag { sample_name }

    input:
    tuple val (sample_name), path (protein_ligand_csv)

    output:
    path ("${sample_name}/index*"), emit: predictions
    path ("${sample_name}/*.npy"), emit: stats

    script:
    """
    #!/bin/bash
    
    shopt -s extglob
    ln -s ${params.diffdock_location}/* .
    shopt -u extglob
    
    mkdir data_local
    
    python datasets/esm_embedding_preparation.py --protein_ligand_csv ${protein_ligand_csv} --out_file data_local/prepared_for_esm.fasta
    HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data_local/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok
    
    python -m inference --protein_ligand_csv ${protein_ligand_csv} --out_dir ${sample_name} \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise
    """
}
