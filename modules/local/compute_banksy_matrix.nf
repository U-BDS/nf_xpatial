process COMPUTE_BANKSY_MATRIX {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(spe_obj), val(k_geom)

    output:
    tuple val(meta), path("*.rds"), val(k_geom), emit: banksy_mtx_spe_obj

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
        --k_geom $k_geom
    """
}
