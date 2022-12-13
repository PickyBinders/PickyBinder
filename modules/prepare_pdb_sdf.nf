/*
*  prepare_pdb_sdf module 
*/

params.OUTPUT = ""

process prepare_pdb_sdf {
    publishDir(params.OUTPUT, mode: 'copy')

    input:
    path (dataset)
    path (output_folder)

    output:
    path "${output_folder}/*.pdb"; emit: pdb_files
    path "${output_folder}/*.sdf"; emit: sdf_files
    path "${output_folder}/protein_ligand.csv"; emit: protein_ligand_csv

    script:
    """
    prepare_sdf_pdf.py ${dataset} ${output_folder}
    """
}
