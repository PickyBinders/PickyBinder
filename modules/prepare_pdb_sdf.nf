/*
*  prepare_pdb_sdf module 
*/

params.OUTPUT = ""

process prepare_pdb_sdf {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("assembly/results/${sample_id}/1_unicycler/last_pilon/pilon", mode: 'copy')
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path(bam), path(assembly)

    output:
    tuple val (sample_id), path("${sample_id}.fasta"), emit: assembly
    tuple val (sample_id), path("${sample_id}.vcf"), emit: vcf
    path "pilon_version.txt", emit: version

    script:
    """
    export _JAVA_OPTIONS="-Xmx10g"
    pilon --mindepth 5 --minmq 10 --threads ${task.cpus} \
    --genome ${assembly} --frags ${bam} --changes --variant \
    --outdir ./ --output ${sample_id}
    pilon --version | cut -d\\  -f1-3 > pilon_vers.txt
    echo ${params.CONTAINER} > pilon_singularity.txt
    cat pilon_vers.txt pilon_singularity.txt | tr "\n" "\t" > pilon_version.txt
    """
}
