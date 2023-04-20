# PickyBinder
Protein-ligand interaction prediction benchmarking workflow

## Overview

This is a nextflow pipeline. Information about Nextflow can be found here https://www.nextflow.io/ 
and in the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html).   

The workflow of the pipeline is defined in the **main.nf** file. The underlying processes 
are defined in the **modules** of each tool. The configuration for the usage of 
the containers (Singularity) and the resource management for the executor (SLURM) 
are defined in the **nextflow.config** file. The locations of the required databases 
and default values are specified in the **params.config** file. 

Each process of the pipeline has its own working directory that is located in 
the **work** folder, where all the output and log files for each process are stored. 
The favoured output files are copied to the results' folder, and therefore it is 
useful to **remove the work directory** if the pipeline finished satisfactorily.  

While the pipeline is running the status can be monitored in **.nextflow.log** or 
in the slurm file. With successful completion the **report.html** is produced which 
gives information about each process including the used resources. 

## Running the workflow

The general command to run the pipeline is:

```
sbatch /scicore/home/schwede/GROUP/TeamLIGATE/tools/PickyBinder/run_pipeline.sh "<Nextflow options> <Workflow options>"
```

The possible options are:

```
Nextflow options:
-with-timeline arg          generates a timeline file at the end of the workflow: <name>.html
-with-dag arg               generates a dag of the workflow: <name>.pdf, <name>.html    

Workflow options:
--pdb_sdf_files	arg         path to pdb and sdf files
--naming arg                naming of the input files: 
                                default (pdbID__ligandName.sdf), 
                                other (ligand needs to have a common identifier with the receptor at 
                                       the start followed by a '_' eg. 6m7h.pdb and 6m7h_ligand.sdf)
--receptor_Hs arg           are the input pdbs hydrated: no (default), yes
--alphafold arg             are the receptors AlphaFold modelled strucures: no (default), yes
--ref_files arg             path to receptor reference files when running workflow with 
                            AlphaFold modelled receptors
--diffdock_mode arg         running DiffDock in batch or single mode: batch (default), single
--autobox_add arg           amount of buffer space to add on each side od the box (default 10)
```


