process RUN_BANKSY {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.5' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.5' 
        }"

    input:
    tuple val(meta), path(xenium_obj)

    output:
    tuple val(meta), path("*.rds"), emit: banksy_xenium_obj

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    run_banksy.R \\
        $args \\
        $assay_flag \\
        --input "$xenium_obj" \\
        --outfile "${prefix}_banksy.rds" \\
        --lambda "${meta.lambda}" \\
        --k_geom "${meta.k_geom}"

    """
}
