process EXTRACT_REDUCED_DIMS {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
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

    def assay_flag = "--assay ${meta.assay}"

    def param_string_flag = meta.clustering_method == 'BANKSY' || meta.clustering_method == 'BANKSYSeurat' ?
        "--param_string l" + "${meta.lambda}" + "-k" + "${meta.k_geom}" + "-d" + "${meta.dim}" :
        "--param_string d" + "${meta.dim}"

    """
    extract_reduced_dims.R \\
        $args \\
        $assay_flag \\
        $param_string_flag \\
        --clustering_method "${meta.clustering_method}" \\
        --input "$xenium_obj" \\
        --outfile "${prefix}_reduced_dims.csv"
    """
}
