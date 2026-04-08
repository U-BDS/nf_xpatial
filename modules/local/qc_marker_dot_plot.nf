process QC_MARKER_DOT_PLOT {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(xenium_obj)
    path marker_list

    output:
    tuple val(meta), path("*.pdf"), emit: dot_plot
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    def cluster_flag = ''
    if ("${meta.clustering_method}" == "BANKSY"){
        cluster_flag = "--cluster_col clust_BSKY_l" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "Seurat"){
        cluster_flag = "--cluster_col clust_SEU_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "BANKSYSeurat"){
        cluster_flag = "--cluster_col clust_BSKYSEU_l" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    }

    """
    qc_marker_dot_plot.R \\
        $args \\
        $cluster_flag \\
        $assay_flag \\
        --input "$xenium_obj" \\
        --marker_list "$marker_list" \\
        --outfile ${prefix}.dot_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
