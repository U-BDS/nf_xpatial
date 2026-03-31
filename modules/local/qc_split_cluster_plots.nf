process QC_SPLIT_CLUSTER_PLOTS {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
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

    def reductions_flag = ''
    if ( "${meta.clustering_method}" == "BANKSY" ){
        reductions_flag = "--reduction BANKSY_umap_l" + "${meta.lambda}" + 
            ".k" + "${meta.k_geom}" + 
            ".d" + "${meta.dim}"
    } else if ( "${meta.clustering_method}" == "Harmony"){
        reductions_flag = "--reduction Harmony_umap_d" + "${meta.dim}"
    } else if ( "${meta.clustering_method}" == "BANKSYSeurat"){
        reductions_flag = "--reduction BANKSYSeurat_umap_d" + "${meta.dim}"
    }

    def cluster_flag = ''
    if ("${meta.clustering_method}" == "BANKSY"){
        cluster_flag = "--cluster_col clust_BSKY_l" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_R" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "Harmony"){
        cluster_flag = "--cluster_col clust_HMY_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "BANKSYSeurat"){
        cluster_flag = "--cluster_col clust_HMY_l" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    }

    """
    qc_split_cluster_plots.R \\
        $args \\
        $reductions_flag \\
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
