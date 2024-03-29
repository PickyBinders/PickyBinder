/*
*  EDMDock module
*/

params.OUTPUT = "$launchDir/predictions/edmdock"

process edmdock {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.edmdock_conda}"

    input:
    path (samples_file)
    path (pdb_Hs)
    path (preped_ligands)
    path (box_files)
    path (edmdock_tool)

    output:
    path("results/**.pdb"), emit: pdbs
    path("results/results.csv"), emit: result_csv

    script:
    """
    mkdir -p dataset local_run
    mkdir -p local_run/paper_baseline
    while read -r line;
    do
        prot=\$(echo \$line | cut -d',' -f1);
	    lig=\$(echo \$line | cut -d',' -f2);
	    box=\$(echo \$line | cut -d',' -f3);
	    mkdir -p dataset/\${box%.box};
	    ln -s \$PWD/\$prot dataset/\${box%.box}/protein.pdb;
	    ln -s \$PWD/\$lig dataset/\${box%.box}/ligand.sdf;
	    cat \$box | cut -d' ' -f3 > tmp.txt;
	    f=\$(cat tmp.txt);
	    echo "\${f//\$'\n'/,}" > dataset/\${box%.box}/box.csv;
    done < edmdock_samples.csv

    cp runs/paper_baseline/config.yml local_run/paper_baseline/config.yml
    cp runs/paper_baseline/weights.ckpt local_run/paper_baseline/weights.ckpt

    export PYTHONPATH="\${PYTHONPATH}:$PWD"

    python scripts/prepare.py --dataset_path dataset
    python scripts/dock.py --run_path local_run/paper_baseline --dataset_path dataset

    mv local_run/paper_baseline/results/ .
    rm -rf local_run
    """
}


process edmdock_single {
    publishDir "$params.OUTPUT/${complex}", mode: 'copy'
    conda "${params.edmdock_conda}"
    tag { "${complex}_${pocket_nr}" }

    input:
    tuple val (complex), val (pocket_nr), val (receptor), val (ligand), path (pdb_Hs), path (preped_ligand), path (box_file)
    path (edmdock_tool)

    output:
    tuple val(complex), val (receptor), val (pocket_nr), path("${pocket_nr}/*.pdb")

    script:
    """
    mkdir -p dataset local_run
    mkdir -p dataset/${complex}_${pocket_nr} local_run/paper_baseline
    mv ${pdb_Hs} dataset/${complex}_${pocket_nr}/protein.pdb
    mv ${preped_ligand} dataset/${complex}_${pocket_nr}/ligand.sdf
    cat ${box_file} | cut -d' ' -f3 > tmp.txt
    f=\$(cat tmp.txt)
    echo "\${f//\$'\n'/,}" > box.csv
    mv box.csv dataset/${complex}_${pocket_nr}/

    cp runs/paper_baseline/config.yml local_run/paper_baseline/config.yml
    cp runs/paper_baseline/weights.ckpt local_run/paper_baseline/weights.ckpt

    export PYTHONPATH="\${PYTHONPATH}:$PWD"

    python scripts/prepare.py --dataset_path dataset
    python scripts/dock.py --run_path local_run/paper_baseline --dataset_path dataset

    mv local_run/paper_baseline/results ./${pocket_nr}
    rm -rf local_run
    """
}
