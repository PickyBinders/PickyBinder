# PickyBinder
This repository contains the implementation of the PickyBinder workflow, a workflow to benchmark 
protein-ligand interaction prediction methods.

The workflow assesses the state-of-the-art protein-ligand complex prediction tools 
Autodock Vina, SMINA, GNINA, DiffDock, and TANKBind. It reports the BiSyRMSD and the 
lDDT-PLI scores calculated with OpenStructure for each predicted ligand pose. For docking tools associated with the 
AutoDock family, the workflow conducts besides blind docking also ligand prediction either
for all binding pockets identified by P2Rank or for a specific user-defined binding site.

## Overview

PickyBinder a nextflow pipeline. Information about Nextflow can be found here https://www.nextflow.io/ 
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

- **Autodock Vina**: Get the Autodock Vina v1.2.5 executable from https://github.com/ccsb-scripps/AutoDock-Vina/releases 
and create a Conda environment for the Python bindings as described in the Autodock Vina manual 
(https://autodock-vina.readthedocs.io/en/latest/installation.html). Either install **meeko** 0.4.0 into the same Conda environment, 
or make an own Conda environment for meeko (https://pypi.org/project/meeko/#2.-rdkit-molecule-from-docking-results).
- **SMINA**: Get the Singularity image for SMINA v2020.12.10 (https://hub.docker.com/r/zengxinzhy/smina, tag:1.0).
- **GNINA**: Get the Singularity image for GNINA v1.0.3 (https://hub.docker.com/r/marcus905/gnina-12/tags, digest: ea3dce32d4a5)
or GNINA v1.0.2 (https://hub.docker.com/r/nmaus/gnina , digest: 7087cbf4dafd).
- **DiffDock**: Create a Conda environment according to the setup guide at https://github.com/gcorso/DiffDock (commit 2c7d438).
- **TANKBind**: Get Singularity image "tankbind_py38" (https://hub.docker.com/r/qizhipei/tankbind_py38, digest: 79a46540b547). 

### Other tools
- **P2Rank**: Get P2Rank v2.4 executable from https://github.com/rdk/p2rank/releases .
- **OpenStructure**: Get Singularity image for OpenStructure 2.5.0 from https://git.scicore.unibas.ch/schwede/openstructure/container_registry/7
- **ADFRsuite**: Build Singularity image from the definition file in the Singularity directory. 

### Preparation of the params.config file

Copy the params.config.in to params.config and adjust the following:
- Provide the paths for the **tool locations**, **conda containers**, and **singularity images**.
- P2Rank needs Java (8 to 18) to run. Eventually, provide a command to **load Java** or leave it blank.
- The **cluster options** section gives the possibility to name the queues, partitions, and other desired SLURM options.
Memory and time are already defined. 

## Input definition

The workflow can be run with multiple complexes at once. For each complex a receptor pdb file and a ligand sdf file
must be provided. As reference for the BiSyRMSD and lDDT-PLI scoring with OpenStructure it is best to use the mmCIF file 
of the receptor. But it is also possible to give a pdb file as the reference. 
If the reference file is a pdb file, then the ligand sdf file is used as the reference ligand.

The input for the pipeline can be allocated in two ways, either prepare a csv file containing
the paths to the relevant files and additional information 
or provide directories for receptor, ligand and reference files.  

### CSV file

The csv file needs to have a header row with the following column names: 
**complex_name,receptor_path,ligand_path_sdf,ligand_path_mol2,reference_path,BS** 

| Column           |                         Content                                | Description                                                                                                                                           |
|:-----------------|:----------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------|
| complex_name     | name used to save predictions                             | if empty the names of the receptor and the ligand will be combined                                                                                    |
| receptor_path    | full path to the receptor pdb file                        |                                                                                                                                                       |
| ligand_path_sdf  | full path to the ligand sdf file                          |                                                                                                                                                       |
| ligand_path_mol2 | full path to ligand mol2 file                             | if the ligand preprocessing fails using the sdf file, ligand preprocessing is retried with the mol2 file. Put a `-` if no mol2 file is available.     |
| reference_path   | full path to the reference file for scoring               | provide either a mmCIF or a pdb file; Put `-` to use the receptor pdb file as the receptor reference and the ligand sdf file as the ligand reference. |
| BS               | x-coordinate_y-coordinate_z-coordinate                    | coordinates of the binding site center; leave empty to use binding pockets predicted by P2Rank                                                        |

Use the workflow option ```--data <input>.csv``` to run the workflow. 

### Input directories

To use the option to give simply a directory with the files as input, the files need to have the following naming pattern:

| Receptor        | Ligand                      | Reference           |
|:----------------|:----------------------------|:--------------------|
| pdbID.pdb       | pdbID__ligandName.sdf       | pdbID.cif/pdb       |
| pdbID_Chain.pdb | pdbID_Chain__ligandName.sdf | pdbID_Chain.cif/pdb |

**Receptor** and **ligand** files can be in the same directory or in separate ones. 

Use the workflow option ```--data <path_to_pdb-sdf-files> ``` for one directory with all input files   
or ```--data <path_to_pdb-files>,<path_to_sdf-files> ``` for two separate directories.

The **reference** files need to be in a different directory. Use the workflow option   
```--ref_files <path_to_reference-files> ```

If the input receptor files should be used as the reference, then the workflow option ```--ref_files```
can be omitted. 

The option to give a specific binding site for the docking tools of the AutoDock family is not available 
when using directories for the input definition. 

## Running the workflow

General command to run the pipeline:

```
nextflow run <path-to-PickyBinder-directory>/main.nf -profile slurm -with-report report.html <Nextflow options> <Workflow options>
```

The workflow execution was built to run with SLURM, but it is also possible to run it locally (remove **-profile slurm**)
or to define a profile for another executor in the nextflow.config file 
(check [Nextflow documentation](https://www.nextflow.io/docs/latest/executor.html) for help).

Available options:

```
Nextflow options:
-----------------
-resume                     resume the pipeline
-c                          give a local configuration file instead of the nextflow.config file
-with-timeline arg          generates a timeline file at the end of the workflow: <name>.html
-with-dag arg               generates a dag of the workflow: <name>.pdf, <name>.html    


Workflow options:
-----------------
--data arg                  <input>.csv, <path_to_pdb-sdf-files>, or <path_to_pdb-files>,<path_to_sdf-files>                           
--naming arg                naming of the input files: default, other
                            default: <pdbID_Chain>.pdb, <pdbID_Chain>__<ligandName>.sdf, <pdbID_Chain>.cif/pdb

--ref_files arg             path to reference files 
                                                 
--receptor_Hs arg           are the input pdbs hydrated: no (default), yes
--alphafold arg             are the receptors AlphaFold modelled strucures: no (default), yes

--tools arg                 comma-separated list of the docking tools to run, default is to 
                            run all available tools:
                            --> diffdock,tankbind,vina,smina,gnina
                            
--scoring_receptors         Compares the receptor structure to the reference using OpenStructure: 
                                no (default), yes (lDDT, RMSD, and QS-score)
--scoring_ligands           Scoring ligand prediction using OpenStructure (lDDT-PLI, BiSyRMSD):
                                yes (default), no                                                   

--diffdock_mode arg         running DiffDock in batch or single mode: batch (default), single
--autobox_add arg           amount of buffer space to add on each side od the box (default 10)
```

## Outputs

The BiSyRMSD and lDDT-PLI of all predicted poses can be found in the ligand_score_summary.csv in the ```scores``` directory. 

All predicted ligand poses can be found in the ```predictions``` directory. 
