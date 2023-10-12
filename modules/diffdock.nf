/*
*  diffdock module
*/

params.OUTPUT = "$launchDir/predictions/diffdock"


process diffdock {
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "diffdock_*.log"
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "protein_ligand_csv_*.csv"
    publishDir "$params.OUTPUT", mode: 'copy', pattern: "diffdock_predictions/**"
    conda "${params.diffdock_conda}"

    input:
    path (protein_ligand_csv)
    path (pdb_files)
    path (sdf_files)
    path (diffd_tool)

    output:
    path ("diffdock_predictions/**"), optional: true, emit: predictions
    path ("diffdock_*.log"), emit: log
    path ("protein_ligand_csv_*.csv"), optional: true

    script:
    """
    time_date=\$(date +"%y-%m-%d-%T")

    if ls ${params.diffdock_tool}/.so3* 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.so3* .; fi
    if ls ${params.diffdock_tool}/.p.npy 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.p.npy .; fi
    if ls ${params.diffdock_tool}/.score.npy 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.score.npy .; fi

    # check which protein-ligand pairs have been docked already and remove them from the protein_ligand.csv
    DIR="$launchDir/predictions/diffdock/diffdock_predictions"
    if [ -d "\$DIR" ]; then
        ls \$DIR > done.txt
        awk -F, '(NR==FNR){a[\$1];next}!(\$1 in a)' done.txt ${protein_ligand_csv} > complexes_to_redo.csv
    else
        cp ${protein_ligand_csv} complexes_to_redo.csv
    fi

    if [[ -s complexes_to_redo.csv ]]; then
        echo 'complex_name,protein_path,ligand_description,protein_sequence' > protein_ligand_withHeader.csv
        cat complexes_to_redo.csv >> protein_ligand_withHeader.csv

        python -m inference --protein_ligand_csv protein_ligand_withHeader.csv --out_dir diffdock_predictions \
            ${params.diffdock_params}

        for dir in diffdock_predictions/*;
        do
            if [ "\$(ls -A \$dir)" ]; then
                for file in \${dir}/*
                do
                    mv \$file \$(dirname \$file)/\$(basename \$dir)_\$(basename \$file)
                done
             fi
        done
        mv protein_ligand_withHeader.csv protein_ligand_csv_\${time_date}.csv
    else
        echo "There are no complexes to dock"
    fi

    ln -s .command.log diffdock_\${time_date}.log
    """
}


process diffdock_single {
    publishDir "$params.OUTPUT/diffdock_predictions", mode: 'copy', pattern: "${complex}/*"
    publishDir "$params.OUTPUT/diffdock_log_files", mode: 'copy', pattern: "diffdock_*.log"
    conda "${params.diffdock_conda}"
    tag { complex }

    input:
    tuple val (ligand), val (receptor), val (complex), path (pdb_file), path (sdf_file)
    path (diffd_tool)

    output:
    tuple val (complex), path ("${complex}/*"), emit: predictions
    path ("diffdock_*.log"), emit: log

    script:
    """
    if ls ${params.diffdock_tool}/.so3* 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.so3* .; fi
    if ls ${params.diffdock_tool}/.p.npy 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.p.npy .; fi
    if ls ${params.diffdock_tool}/.score.npy 1> /dev/null 2>&1; then ln -s ${params.diffdock_tool}/.score.npy .; fi


    python -m inference --complex_name ${complex} --protein_path ${pdb_file} --ligand ${sdf_file} --out_dir ./ \
       ${params.diffdock_params}

    if [ "\$(ls -A ${complex})" ]; then
        for file in ${complex}/*
        do
            mv \$file ${complex}/${complex}_\$(basename \$file)
        done
    fi

    ln -s .command.log diffdock_${complex}.log
    """
}
