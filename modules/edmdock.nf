/*
*  EDMDock module
*/

params.OUTPUT = "$launchDir/predictions/edmdock"

process edmdock {
    publishDir "$params.OUTPUT", mode: 'copy'
    conda "${params.edmdock_conda}"
    tag { complex }

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
    mkdir -p dataset
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

    cp -r runs/ local_run/
    rm local_run/paper_baseline/results/*

    export PYTHONPATH="\${PYTHONPATH}:$PWD"

    python scripts/prepare.py --dataset_path dataset
    python scripts/dock.py --run_path local_run/paper_baseline --dataset_path dataset

    mv local_run/paper_baseline/results/ .
    """
}


process edmdock_single {
    publishDir "$params.OUTPUT/${complex}", mode: 'copy'
    conda "${params.edmdock_conda}"
    tag { complex }

    input:
    tuple val (complex), val (pocket), val (receptor), val (ligand), path (pdb_Hs), path (preped_ligand), path (box_file)
    path (edmdock_tool)

    output:
    path("${pocket}")

    script:
    """
    mkdir dataset
    mkdir dataset/${complex}_${pocket}
    mv ${pdb_Hs} dataset/${complex}_${pocket}/protein.pdb
    mv ${preped_ligand} dataset/${complex}_${pocket}/ligand.sdf
    cat ${box_file} | cut -d' ' -f3 > tmp.txt
    f=\$(cat tmp.txt)
    echo "\${f//\$'\n'/,}" > box.csv
    mv box.csv dataset/${complex}_${pocket}/

    cp -r runs/ local_run/
    rm local_run/paper_baseline/results/*

    export PYTHONPATH="\${PYTHONPATH}:$PWD"

    python scripts/prepare.py --dataset_path dataset
    python scripts/dock.py --run_path local_run/paper_baseline --dataset_path dataset

    mv local_run/paper_baseline/results ./${pocket}
    rm -rf local_run
    """
}
