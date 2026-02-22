process QC_DIM_PLOT_COUNTOUR {
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
    tuple val(meta), path("*.png"), emit: countour_dim_plot
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    def embeddings_flag = ''
    if ( "${meta.clustering_method}" == "BANKSY" ){
        embeddings_flag = "--embedding BANKSY_UMAPBANKSYharmony_d" + "${meta.dim}"
    } else if ( "${meta.clustering_method}" == "Harmony"){
        embeddings_flag = "--embedding Harmony_umap_d" + "${meta.dim}"
    } else if ( "${meta.clustering_method}" == "BANKSYSeurat"){
        embeddings_flag = "--embedding BANKSYSeurat_umap_d" + "${meta.dim}"
    }

    def cluster_flag = ''
    if ("${meta.clustering_method}" == "BANKSY"){
        cluster_flag = "--cluster_col clust_BSKY_AGF1_L" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_R" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "Harmony"){
        cluster_flag = "--cluster_col clust_HMY_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    } else if ("${meta.clustering_method}" == "BANKSYSeurat"){
        cluster_flag = "--cluster_col clust_l" + "${meta.lambda}" + 
            "_k" + "${meta.k_geom}" + 
            "_d" + "${meta.dim}" + 
            "_r" + "${meta.res}"
    }

    """
    qc_dim_plot_countour.R \\
        $args \\
        $embeddings_flag \\
        $cluster_flag \\
        --input "$xenium_obj" \\
        --outfile ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
