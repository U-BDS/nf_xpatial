process EXTRACT_SEURAT_REDUCED_DIMS {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(xenium_obj)

    output:
    tuple val(meta), path("*embeddings*.csv") , emit: embeddings_csv
    tuple val(meta), path("*loadings*.csv")   , emit: loadings_csv, optional: true
    tuple val(meta), path("*stdev*.csv")      , emit: stdev_csv,    optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    def param_string_flag = meta.clustering_method == 'BANKSY' ?
        "--param_string 'l${meta.lambda}_k${meta.k_geom}_n${meta.nPCs}_r${meta.res}'" :
        "--param_string 'd${meta.dim}_r${meta.res}'"

    """
    extract_seurat_reduced_dims.R \\
        $args \\
        $assay_flag \\
        $param_string_flag \\
        --clustering_method "${meta.clustering_method}" \\
        --input "$xenium_obj" \\
        --outfile "${prefix}_reduced_dims.csv"
    """
}
