process CONVERT_SPE_TO_SEURAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(xenium_obj), path(spe_obj)

    output:
    tuple val(meta), path("*.rds"), emit: converted_seurat_object

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    convert_spe_to_seurat.R \\
        $args \\
        $assay_flag \\
        --xenium "$xenium_obj" \\
        --spe_obj "$spe_obj" \\
        --outfile "${prefix}_xenium.rds"
    """
}
