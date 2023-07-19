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
    combine_all_scores.py
    """

}