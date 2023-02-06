/*
*  diffdock module
*/

params.OUTPUT = "$launchDir/diffdock"


process create_diffdock_csv {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'
    
    input: 
    path (ref_sdf_files)
    
    output:
    path("protein_ligand.csv"), emit: csv_for_diffdock
    
    script:
    """
    create_diffdock_csv.py ${ref_sdf_files}
    """
    
}


process diffdock {
    publishDir(params.OUTPUT, mode: 'copy', saveAs: { filename -> if (filename == ".command.log") "diffdock.log"})
    publishDir(params.OUTPUT, mode: 'copy', pattern: "diffdock_predictions") 
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/diffdock'
    label 'diffdock'

    input:
    path (protein_ligand_csv)
    path (pdb_files)
    path (sdf_files)
    path (diffd_tool)

    output:
    path ("diffdock_predictions/"), emit: predictions
    path (".command.log"), emit: diffdock_log

    script:
    """
    mkdir data_local
    
    python datasets/esm_embedding_preparation.py --protein_ligand_csv ${protein_ligand_csv} --out_file data_local/prepared_for_esm.fasta
    HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data_local/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok
    
    python -m inference --protein_ligand_csv ${protein_ligand_csv} --out_dir diffdock_predictions \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise --keep_local_structures
    """
}


process diffdock_single {
    //publishDir(params.OUTPUT, mode: 'copy')
    publishDir("diffdock_predictions", mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/diffdock'
    label 'diffdock'
    tag { complex }

    input:
    tuple val (ligand), val (receptor_chain), val (complex), path (pdb_file), path (sdf_file)
    path (diffd_tool)

    output:
    path ("${complex}/index*"), emit: predictions
    path ("${complex}/*.npy"), emit: stats

    script:
    """
    mkdir data_local
    
    python datasets/esm_embedding_preparation.py --protein_path ${pdb_file} --out_file data_local/prepared_for_esm.fasta
    HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data_local/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok
    
    python -m inference --protein_path ${pdb_file} --ligand ${sdf_file} --out_dir ${complex} \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise --keep_local_structures
    """
}
