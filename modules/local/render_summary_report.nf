process RENDER_SUMMARY_REPORT {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta),
        path(report_template),
        path(filtered_stats_file),
        path(filtered_dim_plot),
        path(filtered_vln_plot),
        path(filtered_feature_scatter_plot),
        path(filtered_nfeature_plot),
        path(filtered_ncount_plot),
        path(cell_shape_proportion_plot),
        path(cell_segmentation_proportion_plot),
        path(compiled_histogram_plot),
        path(compiled_filtered_box_plot),
        path(compiled_filtered_overlapping_histogram_plot),
        path(log_norm_seurat_umap_video),
        path(log_norm_seurat_split_cluster_video),
        path(log_norm_seurat_vln_video),
        path(log_norm_seurat_dot_video),
        path(area_norm_seurat_umap_video),
        path(area_norm_seurat_split_cluster_video),
        path(area_norm_seurat_vln_video),
        path(area_norm_seurat_dot_video),
        path(log_norm_banksy_umap_video),
        path(log_norm_banksy_split_cluster_video),
        path(log_norm_banksy_vln_video),
        path(log_norm_banksy_dot_video),
        path(area_norm_banksy_umap_video),
        path(area_norm_banksy_split_cluster_video),
        path(area_norm_banksy_vln_video),
        path(area_norm_banksy_dot_video)

    output:
    tuple val(meta), path("*.rds"), emit: all_cluster_xenium_obj
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """

    library(knitr)

    rmarkdown::render(
        "${report_template}",
        params = list(
            filtered_stats_file = "${filtered_stats_file.join(',')}",
            filtered_dim_plot = "${filtered_dim_plot}",
            filtered_vln_plot = "${filtered_vln_plot}",
            filtered_feature_scatter_plot = "${filtered_feature_scatter_plot}",
            filtered_nfeature_plot = "${filtered_nfeature_plot}",
            filtered_ncount_plot = "${filtered_ncount_plot}",
            cell_shape_proportion_plot = "${cell_shape_proportion_plot}",
            cell_segmentation_proportion_plot = "${cell_segmentation_proportion_plot}",
            compiled_histogram_plot = "${compiled_histogram_plot}",
            compiled_filtered_box_plot = "${compiled_filtered_box_plot}",
            compiled_filtered_overlapping_histogram_plot = "${compiled_filtered_overlapping_histogram_plot}",
            log_norm_seurat_umap_video = "${log_norm_seurat_umap_video}",
            log_norm_seurat_split_cluster_video = "${log_norm_seurat_split_cluster_video}",
            log_norm_seurat_vln_video = "${log_norm_seurat_vln_video}",
            log_norm_seurat_dot_video = "${log_norm_seurat_dot_video}",
            area_norm_seurat_umap_video = "${area_norm_seurat_umap_video}",
            area_norm_seurat_split_cluster_video = "${area_norm_seurat_split_cluster_video}",
            area_norm_seurat_vln_video = "${area_norm_seurat_vln_video}",
            area_norm_seurat_dot_video = "${area_norm_seurat_dot_video}",
            log_norm_banksy_umap_video = "${log_norm_banksy_umap_video}",
            log_norm_banksy_split_cluster_video = "${log_norm_banksy_split_cluster_video}",
            log_norm_banksy_vln_video = "${log_norm_banksy_vln_video}",
            log_norm_banksy_dot_video = "${log_norm_banksy_dot_video}",
            area_norm_banksy_umap_video = "${area_norm_banksy_umap_video}",
            area_norm_banksy_split_cluster_video = "${area_norm_banksy_split_cluster_video}",
            area_norm_banksy_vln_video = "${area_norm_banksy_vln_video}",
            area_norm_banksy_dot_video = "${area_norm_banksy_dot_video}"
        )
    )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
