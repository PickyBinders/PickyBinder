/*
*  add_Hs_to_receptor module 
*/

params.OUTPUT = "$launchDir/preprocessing/receptor_Hs"

process add_Hs_to_receptor {
    publishDir "$params.OUTPUT", mode: 'copy'
    //container "${params.adfr_sing}"
    conda "${params.vina_conda}"
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_file)

    output:
    tuple val (receptor), path ("*_Hs.pdb"), emit: pdb_Hs

    script:
    """
    ln ${pdb_file} ${receptor}_input.pdb
    run_reduce.py ${receptor}_input.pdb ${receptor} ${params.adfr_sing} $baseDir/reduce_wwPDB_het_dict.txt
    rm ${receptor}_input.pdb
    """
}

//reduce -DB $baseDir/reduce_wwPDB_het_dict.txt ${pdb_file} > ${receptor}_Hs.pdb
