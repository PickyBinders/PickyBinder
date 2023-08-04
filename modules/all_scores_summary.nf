/*
*  all scores summary module
*/

params.OUTPUT = "$launchDir/scores"

process combine_all_scores {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.meeko_conda}"

    input:
    path (all_files_with_tool_scores)

    output:
    path ("all_scores_summary.csv"), emit: all_score_summary

    script:
    """
    if [[ ! -f ligand_score_summary.csv ]]
    then
        echo 'Tool,Complex,Pocket,Rank,lddt_pli,rmsd,Reference_Ligand,center_x,center_y,center_z' > ligand_score_summary.csv
    fi

    combine_all_scores.py $launchDir
    """

}