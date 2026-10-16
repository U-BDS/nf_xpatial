process QC_SPLIT_CLUSTER_PLOTS {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.5' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.5' 
        }"

    input:
    tuple val(meta), path(xenium_obj)

    output:
    tuple val(meta), path("*.png"), emit: split_cluster_plot
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def assay_flag = "--assay ${meta.assay}"
    def embeddings_flag = "--reduction ${meta.embedding_name}"
    def cluster_flag = "--cluster_col ${meta.cluster_name}"

    """
    qc_split_cluster_plots.R \\
        $args \\
        $embeddings_flag \\
        $cluster_flag \\
        $assay_flag \\
        --input "$xenium_obj" \\
        --outfile ${prefix}_split_cluster_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
