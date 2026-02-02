process RUN_HARMONY_BANKSY {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(spe_obj)

    output:
    tuple val(meta), path("*.rds"), emit: banksy_pca_harmony_obj

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    run_harmony_banksy.R \\
        $args \\
        $assay_flag \\
        --ncores ${task.cpus} \\
        --input "$spe_obj" \\
        --outfile "${prefix}_banksy_harmony_spe.rds"
    """
}
