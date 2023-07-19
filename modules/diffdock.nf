/*
*  diffdock module
*/

params.OUTPUT = "$launchDir/predictions/diffdock"


process diffdock {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "diffdock_*.log"
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "protein_ligand_withHeader.csv"
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "diffdock_predictions/**"
    conda "${params.diffdock_conda}"

    input:
    path (protein_ligand_csv)
    path (pdb_files)
    path (sdf_files)
    path (diffd_tool)

    output:
    path ("diffdock_predictions/**"), emit: predictions
    path ("diffdock_*.log")
    path ("protein_ligand_withHeader.csv")

    script:
    """
    echo 'complex_name,protein_path,ligand_description,protein_sequence' > protein_ligand_withHeader.csv
    cat ${protein_ligand_csv} >> protein_ligand_withHeader.csv

    python -m inference --protein_ligand_csv protein_ligand_withHeader.csv --out_dir diffdock_predictions \
       --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise

    for dir in diffdock_predictions/*; do for file in \${dir}/*; do mv \$file \$(dirname \$file)/\$(basename \$dir)_\$(basename \$file);done;done

    time_date=\$(date +"%y-%m-%d-%T")
    ln -s .command.log diffdock_\${time_date}.log
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

    for file in ${complex}/*; do mv \$file ${complex}/${complex}_\$(basename \$file);done
    """
}
