process EXTRACT_CLUSTER_METADATA {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
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

    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    def param_string_flag = meta.clustering_method == 'BANKSY' ?
        "--param_string l" + "${meta.lambda}" + "_k" + "${meta.k_geom}" + "_n" + "${meta.nPCs}" + "_r" + "${meta.res}" :
        "--param_string d" + "${meta.dim}" + "_r" + "${meta.res}"

    """
    extract_cluster_metadata.R \\
        $args \\
        $assay_flag \\
        $param_string_flag \\
        --input "$spe_obj" \\
        --outfile "${prefix}_clusts.tsv"
    """
}
