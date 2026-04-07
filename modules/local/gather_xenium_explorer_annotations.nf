process GATHER_XENIUM_EXPLORER_ANNOTATIONS {
    tag "$meta.id"
    label 'process_low'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(xenium_explorer_annotations)

    output:
    tuple val(meta), path("*.tsv"), emit: manual_annotations
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def input_file_arg = xenium_explorer_annotations.split(',').collect{ file -> "-i ${file.trim()}"}.join(' ')

    """
    gather_xenium_explorer_annotations.sh \\
        $args \\
        $input_file_arg \\
        -o ${prefix}_manual_annotations.tsv
    """
}
