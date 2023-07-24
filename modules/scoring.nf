/*
*  scoring module 
*/

params.OUTPUT = "$launchDir/scores"

process pdb_to_sdf {
    publishDir "${out_dir}", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"

    input:
    path(prediction_pdb)
    val(out_dir)

    output:
    path ("*.sdf")

    script:
    """
    python3 $baseDir/bin/pdb_to_sdf.py ${prediction_pdb}
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
            ost compare-ligand-structures --substructure-match -m ${model_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} -o \${model%.sdf}.json --lddt-pli --rmsd
        elif [[ ${ref_receptor} == *.cif ]]
        then
            ost compare-ligand-structures --substructure-match -m ${model_receptor} -ml \${model} -r ${ref_receptor} -o \${model%.sdf}.json --lddt-pli --rmsd
        fi
    done

    python3 $baseDir/bin/combine_ost_scores.py ${tool_name} ${complex}
    """
}


process score_summary {
    publishDir "$params.OUTPUT", mode: 'copy'

    input:
    path (scores)

    output:
    path ("ligand_score_summary.csv"), emit: score_summary

    script:
    """
    if [[ ! -f $launchDir/scores/ligand_score_summary.csv ]]
    then
        echo 'Tool,Complex,Pocket,Rank,lddt_pli,rmsd,Reference_Ligand,center_x,center_y,center_z' > score_summary.csv
        grep -v 'Tool' *_score_summary.csv | cut -d':' -f2 >> score_summary.csv
    else
        cp $launchDir/scores/ligand_score_summary.csv score_summary.csv
        grep -v 'Tool' *_score_summary.csv | cut -d':' -f2 >> score_summary.csv
    fi

    (head -n 1 score_summary.csv && tail -n +2 score_summary.csv | sort) | uniq > ligand_score_summary.csv
    """
}


process ost_scoring_receptors {
    publishDir "$params.OUTPUT/receptors", mode: 'copy'
    container "${params.ost_sing}"
    containerOptions "-B $baseDir/bin"
    tag { receptor }

    input:
    tuple val (receptor), path (ref_receptor), path (model_receptor)

    output:
    tuple val (receptor), path ("*.json"), emit: score

    script:
    """
    ost compare-structures -m ${model_receptor} -r ${ref_receptor} -o ${receptor}.json --lddt --qs-score --rigid-scores
    """
}


process combine_receptors_scores {
    publishDir "$params.OUTPUT", mode: 'copy'

    input:
    path (scores)

    output:
    path ("*.csv"), emit: score_summary

    script:
    """
    echo receptor, lddt, rmsd, qs_global > receptor_score_summary.csv

    for file in *.json
    do
        lddt=\$(grep 'lddt' \$file | cut -d' ' -f6)
        rmsd=\$(grep 'rmsd' \$file | cut -d' ' -f6)
        qs_global=\$(grep 'qs_global' \$file | cut -d' ' -f6)
        echo \${file%.json}, \$lddt \$rmsd \$qs_global >> receptor_score_summary.csv
    done
    """
}
