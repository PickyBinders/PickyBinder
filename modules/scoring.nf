/*
*  scoring module 
*/

params.CONTAINER = "registry.scicore.unibas.ch-schwede-openstructure-2.4.0"
params.OUTPUT = "$launchDir/scores"

process ost_scoring {
    publishDir "$params.OUTPUT/${complex}/${tool_name}", mode: 'copy'
    container params.CONTAINER
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor), path (ref_ligand), path (modeled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score
    tuple val (complex), path ("*.csv"), emit: summary

    script:
    """
    for model in ${modeled_ligands}
        do
        ost compare-ligand-structures -m ${ref_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} -o \${model%.sdf}.json --lddt-pli --rmsd
        done

    python3 $baseDir/bin/combine_ost_scores.py ${tool_name}
    """
}


process ost_scoring_modelLigands {
    publishDir "$params.OUTPUT/${complex}/${tool_name}", mode: 'copy'
    container params.CONTAINER
    containerOptions "-B $baseDir/bin"
    tag { complex }

    input:
    tuple val (complex), val (receptor), val (ligand), path (ref_receptor), path (model_receptor), path (ref_ligand), path (modeled_ligands)
    val (tool_name)

    output:
    tuple val (complex), path ("*.json"), emit: score
    tuple val (complex), path ("*.csv"), emit: summary

    script:
    """
    for model in ${modeled_ligands}
        do
        ost compare-ligand-structures -m ${model_receptor} -ml \${model} -r ${ref_receptor} -rl ${ref_ligand} -o \${model%.sdf}.json --lddt-pli --rmsd
        done

    python3 $baseDir/bin/combine_ost_scores.py ${tool_name}
    """
}


process ost_scoring_modelReceptors {
    publishDir "$params.OUTPUT/receptors", mode: 'copy'
    container params.CONTAINER
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


process combine_modelReceptors_scores {
    publishDir "$params.OUTPUT/receptors", mode: 'copy'

    input:
    path (scores)

    output:
    path ("*.csv"), emit: score_summary

    script:
    """
    echo receptor, lddt, rmsd, qs_global > modelled_receptors_score_summary.csv

    for file in *.json
    do
        lddt=\$(grep 'lddt' \$file | cut -d' ' -f6)
        rmsd=\$(grep 'rmsd' \$file | cut -d' ' -f6)
        qs_global=\$(grep 'qs_global' \$file | cut -d' ' -f6)
        echo \${file%.json}, \$lddt \$rmsd \$qs_global >> modelled_receptors_score_summary.csv
    done
    """
}

process score_summary {
    publishDir "$params.OUTPUT", mode: 'copy'

    input:
    path (scores)

    output:
    path ("score_summary.csv"), emit: score_summary

    script:
    """
    echo 'Tool,Complex,Pocket,Rank,lddt_pli,rmsd' > score_summary.csv
    grep -v 'Tool' *_score_summary.csv | cut -d':' -f2 >> score_summary.csv
    """
}
