/*
*  scoring module 
*/

params.OUTPUT = "$launchDir/scores"

process pdb_to_sdf_batch {
    publishDir "${out_dir}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"

    input:
    path(prediction_pdb)
    val(out_dir)

    output:
    path ("*.sdf"), emit: sdf_files
    path ("pdb_to_sdf*.log")

    script:
    """
    python3 $baseDir/bin/pdb_to_sdf.py ${prediction_pdb}

    run_name=${params.timestamp}
    ln -s .command.log pdb_to_sdf_\${run_name}.log
    """
}


process pdb_to_sdf_single {
    publishDir "${out_dir}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val(complex), val (receptor), path(prediction_pdb)
    val(out_dir)

    output:
    path ("*.sdf"), emit: sdf_files

    script:
    """
    python3 $baseDir/bin/pdb_to_sdf.py ${prediction_pdb}
    """
}


process ost_scoring_single {
    publishDir "$params.OUTPUT/ligands/", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor, stageAs: "ref/*"), path (model_receptor), path (model_ligand)

    output:
    tuple val (complex), path ("*.json"), emit: score

    script:
    """
    ost compare-ligand-structures --substructure-match \
        -m ${model_receptor} -ml ${model_ligand} -r ${ref_receptor} \
        -o ${complex}.json \
        --lddt-pli --rmsd
    """
}


process ost_scoring_single_summary {
    publishDir "$params.OUTPUT", mode: 'copy'
    //container "${params.ost_sing}"
    //containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    path (scores)

    output:
    path ("*.csv"), emit: score_summary

    script:
    """
    source /scicore/home/schwede/leeman0000/activate-ost-develop

    python3 $baseDir/bin/single_ost_scores_summary.py
    """
}



process ost_scoring {
    publishDir "$params.OUTPUT/ligands/${complex}/${tool_name}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor, stageAs: "ref/*"), path (model_receptor), path (ref_ligand), path (modelled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score
    tuple val (complex), val (receptor), path ("*.csv"), emit: summary

    script:
    """
    for model in ${modelled_ligands}
    do
        if [[ ${ref_receptor} == *.pdb ]]
        then
            ost compare-ligand-structures --substructure-match \
                -m ${model_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} \
                -o \${model%.sdf}.json \
                --lddt-pli --rmsd \
                || continue
        elif [[ ${ref_receptor} == *.cif ]]
        then
            ost compare-ligand-structures --substructure-match \
                -m ${model_receptor} -ml \${model} -r ${ref_receptor} \
                -o \${model%.sdf}.json \
                --lddt-pli --rmsd \
                || continue
        fi
    done

    python3 $baseDir/bin/combine_ost_scores.py ${tool_name} ${complex}
    """
}


process ost_scoring_diffdock {
    publishDir "$params.OUTPUT/ligands/${complex}/${tool_name}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor, stageAs: "ref/*"), path (model_receptor), path (ref_ligand), path (modelled_ligands)
    val (tool_name)

    output:
    tuple val (complex), val(receptor), path ("*.json"), emit: scores

    script:
    """
    for model in ${modelled_ligands}
    do
        if [[ ${ref_receptor} == *.pdb ]]
        then
            ost compare-ligand-structures --substructure-match \
                -m ${model_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} \
                -o \${model%.sdf}.json \
                --lddt-pli --rmsd \
                || continue
        elif [[ ${ref_receptor} == *.cif ]]
        then
            ost compare-ligand-structures --substructure-match \
                -m ${model_receptor} -ml \${model} -r ${ref_receptor} \
                -o \${model%.sdf}.json \
                --lddt-pli --rmsd \
                || continue
        fi
    done
    """
}


process combine_dd_scores {
    publishDir "$params.OUTPUT/ligands/${complex}/${tool_name}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val(receptor), path (score_files)
    val (tool_name)

    output:
    tuple val (complex), val (receptor), path ("*.csv"), emit: summary

    script:
    """
    python3 $baseDir/bin/combine_ost_scores.py ${tool_name} ${complex}
    """
}


process ost_score_summary {
    publishDir "$params.OUTPUT/summary_files_$params.runID", mode: 'copy'

    input:
    path (scores)

    output:
    path ("ligand_score_summary.csv"), emit: score_summary

    script:
    """
    echo 'Tool,Complex,Pocket,Rank,lDDT-PLI,lDDT-LP,BiSyRMSD,Reference_Ligand,Box_Center_x,Box_Center_y,Box_Center_z' > summary.csv
    grep -v 'Tool' *_score_summary.csv | cut -d':' -f2 >> summary.csv

    mv summary.csv ligand_score_summary.csv
    """
}


process ost_scoring_receptors {
    publishDir "$params.OUTPUT/receptors", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { receptor }

    input:
    tuple val (receptor), path (ref_receptor, stageAs: "ref/*"), path (model_receptor)

    output:
    tuple val (receptor), path ("*.json"), emit: score

    script:
    """
    if [[ ${ref_receptor} == *.pdb ]]
    then
        ost compare-structures -m ${model_receptor} -r ${ref_receptor} -o ${receptor}_pdbRef.json --lddt --qs-score --rigid-scores
    elif [[ ${ref_receptor} == *.cif ]]
    then
        ost compare-structures -m ${model_receptor} -r ${ref_receptor} -o ${receptor}.json --lddt --qs-score --rigid-scores
    fi
    """
}


process combine_receptors_scores {
    publishDir "$params.OUTPUT/summary_files_$params.runID", mode: 'copy'

    input:
    path (scores)

    output:
    path ("*.csv"), emit: score_summary

    script:
    """
    echo receptor, lddt, rmsd, qs_global > receptor_score_summary.csv

    for file in *.json
    do
        lddt=\$(grep '"lddt"' \$file | cut -d' ' -f6)
        rmsd=\$(grep '"rmsd"' \$file | cut -d' ' -f6)
        qs_global=\$(grep '"qs_global"' \$file | cut -d' ' -f6)
        echo \${file%.json}, \$lddt \$rmsd \$qs_global >> receptor_score_summary.csv
    done
    """
}
