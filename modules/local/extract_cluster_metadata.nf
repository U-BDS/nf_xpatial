process EXTRACT_CLUSTER_METADATA {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(spe_obj)

    output:
    tuple val(meta), path("*.tsv"), emit: cluster_metadata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def assay_flag = "--assay ${meta.assay}"

    def param_string_flag = meta.clustering_method == 'BANKSY' || meta.clustering_method == 'BANKSYSeurat' ?
        "--param_string l" + "${meta.lambda}" + "_k" + "${meta.k_geom}" + "_d" + "${meta.dim}" + "_r" + "${meta.res}" :
        "--param_string d" + "${meta.dim}" + "_r" + "${meta.res}"

    """
    extract_cluster_metadata.R \\
        $args \\
        $assay_flag \\
        $param_string_flag \\
        --clustering_method "${meta.clustering_method}" \\
        --input "$spe_obj" \\
        --outfile "${prefix}_clusts.tsv"
    """
}
