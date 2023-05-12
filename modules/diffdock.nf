/*
*  diffdock module
*/

params.OUTPUT = "$launchDir/diffdock"

process create_diffdock_csv {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.meeko_conda}"

    input:
    path (ref_sdf_files)

    output:
    path("protein_ligand.csv"), emit: csv_for_diffdock

    script:
    """
    create_diffdock_csv_new.py ${params.naming} ${params.receptor_Hs} ${ref_sdf_files}
    """
}


process diffdock {
    publishDir "$params.OUTPUT", mode: 'copy', saveAs: { filename -> if (filename == ".command.log") "diffdock.log"}
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "diffdock_predictions/**"
    conda "${params.diffdock_conda}"

    input:
    path (protein_ligand_csv)
    path (pdb_files)
    path (sdf_files)
    path (diffd_tool)

    output:
    path ("diffdock_predictions/**"), emit: predictions
    path (".command.log"), emit: diffdock_log

    script:
    """
    if [ ${params.naming} != default ]
    then
        for file in ${pdb_files}; do receptor=\$(basename \$file .pdb | cut -d'_' -f1); mv \$file \${receptor}_receptor.pdb;done
    fi

    python -m inference --protein_ligand_csv ${protein_ligand_csv} --out_dir diffdock_predictions \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise
    """
}


process diffdock_single {
    publishDir "$params.OUTPUT/diffdock_predictions", mode: 'copy'
    conda "${params.diffdock_conda}"
    tag { complex }

    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_file), path (sdf_file)
    path (diffd_tool)

    output:
    tuple val (complex), path ("${complex}/*"), emit: predictions

    script:
    """
    python -m inference --complex_name ${complex} --protein_path ${pdb_file} --ligand ${sdf_file} --out_dir ./ \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise
    """
}
