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
    path ("*_summary.csv"), optional: true

    script:
    """
    if [[ ! -f ligand_score_summary.csv ]]
    then
        echo 'Tool,Complex,Pocket,Rank,lDDT-PLI,lDDT-LP,BiSyRMSD,Reference_Ligand,Box_Center_x,Box_Center_y,Box_Center_z' > ligand_score_summary.csv
    fi

    combine_all_scores.py $launchDir
    """

}