process EXTRACT_SEURAT_CLUSTER_METADATA {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(spe_obj), val(dim), val(res)

    output:
    tuple val(meta), path("*.tsv"), emit: cluster_metadata

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    extract_seurat_cluster_metadata.R \\
        $args \\
        $assay_flag \\
        --input "$spe_obj" \\
        --dim $dim \\
        --res $res \\
        --outfile "${prefix}_clusts.tsv"
    """
}
