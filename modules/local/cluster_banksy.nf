process CLUSTER_BANKSY {
    tag "$meta.id"
    label 'process_low'

    //container "nf_xenium_analysis_0.0.1.sif"
    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.1' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.1' 
        }"

    input:
    tuple val(meta), path(spe_obj), val(nPCs), val(lambda), val(res)

    output:
    tuple val(meta), path("*.rds"), emit: banksy_cluster_spe_obj

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    cluster_banksy.R \\
        $args \\
        $assay_flag \\
        --input "$spe_obj" \\
        --outfile "${prefix}_banksy_pca_spe.rds" \\
        --lambda $lambda \\
        --nPCs $nPCs \\
        --res $res

    """
}
