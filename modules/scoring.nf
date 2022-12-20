/*
*  scoring module 
*/

params.OUTPUT = "$launchDir"

process rmsd {
    publishDir(params.OUTPUT, mode: 'copy')
    conda '/scicore/home/schwede/leeman0000/miniconda3/envs/spyrmsd'

    input:
    path (diffdock_predictions)

    output:
    path ("${params.rmsd_output}/*_performance.png"), emit: performance_plot
    path ("${params.rmsd_output}/*_RMSD.png"), emit: RMSD_plot
    path ("${params.rmsd_output}/*_RMSD.csv"), emit: RMSD_table
    path ("${params.rmsd_output}/summary_table.csv"), emit: summary_table

    script:
    """
    mkdir -p ${params.rmsd_output}
    if [[ ! -d diffdock_predictions ]]; then mkdir diffdock_predictions && mv ${diffdock_predictions} diffdock_predictions; fi
    
    rmsd_scoring.py ${params.pdb_sdf_files} ./diffdock_predictions ${params.rmsd_output}
    """
}
