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
    path ("${params.rmsd_output}/*.csv"), emit: RMSD_table

    script:
    """
    mkdir -p ${params.rmsd_output}
    
    rmsd_scoring.py ${params.pdb_sdf_files} $PWD/diffdock_predictions ${params.rmsd_output}
    """
}
