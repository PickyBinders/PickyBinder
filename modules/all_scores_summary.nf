/*
*  all scores summary module
*/

params.OUTPUT = "$launchDir/scores"

process combine_all_scores {
    publishDir "$params.OUTPUT/summary_files_$params.runID", mode: 'copy'
    conda "${params.meeko_conda}"

    input:
    path (all_files_with_tool_scores)
    path (coordinates)

    output:
    path ("*_summary.csv"), optional: true
    val true, emit: ready

    script:
    """
    if [[ ! -f ligand_score_summary.csv ]]
    then
        echo 'Tool,Complex,Pocket,Rank,lDDT-PLI,lDDT-LP,BiSyRMSD,Reference_Ligand,Box_Center_x,Box_Center_y,Box_Center_z' > ligand_score_summary.csv
    fi

    combine_all_scores.py $launchDir

    ligand_scoring=${params.scoring_ligands}
    if [[ \${ligand_scoring,,} == no ]]; then rm ligand_score_summary.csv; fi
    """

}
