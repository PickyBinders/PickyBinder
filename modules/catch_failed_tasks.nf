/*
*  error reports module
*/

params.OUTPUT = "$launchDir/errors_and_problems"

process catch_ignored_tasks {

    input:
    val (ready)

    output:
    path ("*.csv"), optional: true, emit: failed_tasks

    script:
    """
    trace_report=\$(ls -t $launchDir/PickyBinder_trace_report_*.csv | head -1)
    echo 'Process,Exit_status,Work_directory' > ignored_tasks_\$(basename \$trace_report)
    grep IGNORE \$trace_report | cut -d',' -f1,4,20 >> ignored_tasks_\$(basename \$trace_report) || continue
    """
}


process catch_diffdock_problems {
    publishDir "$params.OUTPUT/$params.runID", mode: 'copy'

    input:
    path (diffdock_log)

    output:
    path ("diffdock*_problems.txt"), optional: true, emit: diffdock_problems

    script:
    """
    log_file=${diffdock_log}
    outfile=\${log_file%.log}_problems.txt

    grep 'Failed for' $diffdock_log > \${outfile} || continue
    grep 'Failed on' $diffdock_log >> \${outfile} || continue

    echo -e "\n" >> \${outfile}

    grep 'Skipped' $diffdock_log >> \${outfile} || continue
    grep -A 1 'Skipping' $diffdock_log | grep -v HAPPENING >> \${outfile} || continue
    grep 'We are skipping it' $diffdock_log >> \${outfile} || continue
    """
}


process no_box_size {
    publishDir "$params.OUTPUT/$params.runID", mode: 'copy'

    input:
    path (box_size_logs)

    output:
    path ("errors_box_size*.txt"), optional: true, emit: errors_box_size

    script:
    """
    run_name=${params.timestamp}

    for file in ${box_size_logs}
    do
        grep -B 3 'The box size is set to' \$file >> errors_box_size.txt || continue
        echo -e "\n" >> errors_box_size.txt
    done

    awk 'NF {p=1} p' errors_box_size.txt | tac | awk 'NF {p=1} p' | tac > errors_box_size_\${run_name}.txt || continue
    rm errors_box_size.txt
    """
}


process error_and_problems_summary {
    publishDir "$params.OUTPUT/$params.runID", mode: 'copy'

    input:
    tuple path (error_problem_files), val (p2rank_done)

    output:
    path ("errors_and_problems_summary_*.txt")

    script:
    """
    run_name=${params.timestamp}

    echo -e "*** Failed processes ***" > errors_and_problems_summary_\${run_name}.txt
    cat ignored_tasks_*.csv >> errors_and_problems_summary_\${run_name}.txt

    if ls ligand_preprocessing_*.log 1> /dev/null 2>&1; then
        echo -e "\n\n*** Ligand preprocessing problems ***" >> errors_and_problems_summary_\${run_name}.txt
        cat ligand_preprocessing_*.log >> errors_and_problems_summary_\${run_name}.txt
    fi

    if ls errors_box_size_*.txt 1> /dev/null 2>&1; then
        echo -e "\n\n*** Errors in box size calculation ***" >> errors_and_problems_summary_\${run_name}.txt
        echo -e "For the following ligands the box size was set to 40Å:" >> errors_and_problems_summary_\${run_name}.txt
        cat errors_box_size_*.txt >> errors_and_problems_summary_\${run_name}.txt
    fi

    if [[ ${p2rank_done} == "done" ]]; then
        echo -e "\n\n*** No binding pocket detected by P2Rank ***" >> errors_and_problems_summary_\${run_name}.txt
        if [[ -f p2rank_no_pockets_found.csv ]]; then
            cat p2rank_no_pockets_found.csv >> errors_and_problems_summary_\${run_name}.txt
        fi
    fi

    if ls diffdock_*_problems.txt 1> /dev/null 2>&1; then
        echo -e "\n\n*** Diffdock failures ***" >> errors_and_problems_summary_\${run_name}.txt
        cat diffdock_*_problems.txt >> errors_and_problems_summary_\${run_name}.txt
    fi

    echo -e "\n" >> errors_and_problems_summary_\${run_name}.txt
    """
}