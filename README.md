# PickyBinder
This repository contains the implementation of the PickyBinder workflow, 
a workflow to benchmark protein-ligand interaction prediction methods.

## Overview

This is a nextflow pipeline. Information about Nextflow can be found here https://www.nextflow.io/ 
and in the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html).   

The workflow of the pipeline is defined in the **main.nf** file. The underlying processes 
are defined in the **modules** of each tool. The configuration for the usage of 
the containers (Singularity) and the resource management for the executor (SLURM) 
are defined in the **nextflow.config** file. The locations of the required tools, Singularity images 
and default values are specified in the **params.config** file. 

Each process of the pipeline has its own working directory that is located in 
the **work** folder, where all the output and log files for each process are stored. 
The designated output files are automatically copied to the results' folder, and therefore it is 
useful to **remove the work directory** if the pipeline finished satisfactorily.  

Upon successful completion the **report.html** is produced which 
gives information about each process including the used resources. 

## Dependencies

### Nextflow

The workflow has been tested using Nextflow 20.10.0. Get Nextflow from https://www.nextflow.io/ .

### Protein-ligand prediction tools

- **Autodock Vina**: Get the Autodock Vina v1.2.3 executable from https://github.com/ccsb-scripps/AutoDock-Vina/releases 
and create a Conda environment for the Python bindings as described in the Autodock Vina manual 
(https://autodock-vina.readthedocs.io/en/latest/installation.html). Either install **meeko** 0.4.0 into the same Conda environment, 
or make an own Conda environment for meeko (https://pypi.org/project/meeko/#2.-rdkit-molecule-from-docking-results).
- **SMINA**: Creat a Conda environment for SMINA v2020.12.10 (https://anaconda.org/conda-forge/smina). 
- **GNINA**: Get the Singularity image for GNINA v1.0.2 (https://hub.docker.com/r/nmaus/gnina , digest: 7087cbf4dafd).
- **DiffDock**: Create a Conda environment according to the setup guide at https://github.com/gcorso/DiffDock (commit 2c7d438).
- **TANKBind**: Get Singularity image "tankbind_py38" (https://hub.docker.com/r/qizhipei/tankbind_py38, digest: 79a46540b547). 

### Other tools
- **P2Rank**: Get P2Rank v2.4 executable from https://github.com/rdk/p2rank/releases .
- **OpenStrucutre**: Get Singularity image for OpenStructure 2.4.0 from https://git.scicore.unibas.ch/schwede/openstructure/container_registry/7
- **ADFRsuite**: Build Singularity image from the definition file in the Singularity directory. 

### Preparation of the params.config file

Copy the params.config.in to params.config and adjust the following:
- Provide the paths for the **tool locations**, **conda containers**, and **singularity images**.
- P2Rank needs Java (8 to 18) to run. Eventually, provide a command to **load Java** or leave it blank.
- The **cluster options** section gives the possibility to name the queues, partitions, and other desired SLURM options.
Memory and time are already defined. 

## Running the workflow

The workflow takes as input a directory containing receptor.pdb and ligand.sdf files. If the analysis is run on 
AlphaFold modeled receptors, a directory with reference pdb files needs to be provided to do the scoring. 

The general command to run the pipeline:

```
nextflow run <path-to-PickyBinder-directory>/main.nf -profile slurm -with-report report.html <Nextflow options> <Workflow options>
```

The workflow execution was built to run with SLURM, but it is also possible to run it locally (remove **-profile slurm**).

Available options:

```
Nextflow options:
-with-timeline arg          generates a timeline file at the end of the workflow: <name>.html
-with-dag arg               generates a dag of the workflow: <name>.pdf, <name>.html    

Workflow options:
--pdb_sdf_files	arg         path to pdb and sdf files
--naming arg                naming of the input files: 
                                default (pdbID_Chain__ligandName.sdf), 
                                other (ligand needs to have a common identifier with the receptor at 
                                       the start followed by a '_' eg. 6m7h.pdb and 6m7h_ligand.sdf)
--receptor_Hs arg           are the input pdbs hydrated: no (default), yes
--alphafold arg             are the receptors AlphaFold modelled strucures: no (default), yes
--ref_files arg             path to receptor reference files when running workflow with 
                            AlphaFold modelled receptors
--diffdock_mode arg         running DiffDock in batch or single mode: batch (default), single
--autobox_add arg           amount of buffer space to add on each side od the box (default 10)
```


