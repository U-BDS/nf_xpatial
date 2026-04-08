process CONVERT_SEURAT_TO_SPE {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(norm_obj)

    output:
    tuple val(meta), path("*.rds"), emit: spe_object

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    convert_seurat_to_spe.R \\
        $args \\
        $assay_flag \\
        --input "$norm_obj" \\
        --outfile "${prefix}_spe.rds"
    """
}
