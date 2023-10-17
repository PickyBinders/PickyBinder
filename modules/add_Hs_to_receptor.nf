/*
*  add_Hs_to_receptor module 
*/

params.OUTPUT = "$launchDir/preprocessing/receptor_Hs"

process add_Hs_to_receptor {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.meeko_conda}"
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_file), val (alphafold)

    output:
    tuple val (receptor), path ("*_Hs.pdb"), val (alphafold), emit: pdb_Hs

    script:
    """
    cp ${pdb_file} ${receptor}_input.pdb
    run_reduce.py ${receptor}_input.pdb ${receptor} ${params.adfr_sing} $baseDir/reduce_wwPDB_het_dict.txt
    rm ${receptor}_input.pdb
    """
}


process fix_pdb {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.pdbfixer_conda}"
    tag { receptor }

    input:
    tuple val (receptor), path (pdb_file, stageAs: "input.pdb")

    output:
    tuple val (receptor), path ("*_Hs.pdb"), emit: pdb_Hs

    script:
    """
    pdbfixer_fix_pdb.py ${pdb_file} ${receptor}
    """
}
