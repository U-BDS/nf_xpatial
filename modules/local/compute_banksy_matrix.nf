process COMPUTE_BANKSY_MATRIX {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.5' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.5' 
        }"

    input:
    tuple val(meta), path(spe_obj)

    output:
    tuple val(meta), path("*.rds"), emit: banksy_mtx_spe_obj

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    compute_banksy_matrix.R \\
        $args \\
        $assay_flag \\
        --input "$spe_obj" \\
        --outfile "${prefix}_banksy_mtx_spe.rds" \\
        --k_geom "${meta.k_geom}"
    """
}
