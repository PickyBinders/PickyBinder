#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Define the pipeline parameters
 */

log.info """
PickyBinder  ~  version ${workflow.manifest.version}
=============================================

Input data             : ${params.data}

=============================================

"""


/*
* define the input channels
*/

// get input csv file
Channel
    .fromPath(params.data)
    .splitCsv(header:true)
    .map{ row-> tuple(file(row.receptor), file(row.ligand), file(row.reference), row.complex_name ) }
    .branch{
        with_complex_name: it[3] != '' && it[3] != '-'
        other: true
    }
    .set{ all_input }

// define a complex name where it is not given
all_input.other.map{ [ it[0], it[1], it[2], it[0].simpleName + '__' + it[1].simpleName ] }
               .concat(all_input.with_complex_name)
               .set{ all_input_defined }

// map pdb and sdf
all_input_defined.map{ [it[0]] }.flatten().set{ pdb_files }
all_input_defined.map{ [it[1]] }.flatten().set{ ref_sdf_files }

// ref files; if none is given assign the input pdb as reference
all_input_defined.branch{
                    with_ref_file: it[2].simpleName != '' && it[2].simpleName != '-'
                    other: true
                    }
                    .set{ ref_definition }

ref_definition.other.map{ [it[3], it[0].simpleName, it[1].simpleName, it[0], it[0], it[1]]}         // complex, receptor, ligand, ref_receptor_file, model_receptor_file, ref_ligand_file
                    .concat(ref_definition.with_ref_file.map{ [it[3], it[0].simpleName, it[1].simpleName, it[2], it[0], it[1]]})
                    .set{ scoring_ref }

/*
* include the modules
*/

include { ost_scoring_single; ost_scoring_single_summary } from "./modules/scoring"

/*
* main workflow
*/

workflow {

    // scoring
    ost_scores = ost_scoring_single( scoring_ref )
    all_scores = ost_scoring_single_summary( ost_scores.map{ it[1] }.collect() )
}

