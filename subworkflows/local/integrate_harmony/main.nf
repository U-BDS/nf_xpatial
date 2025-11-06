#!/usr/bin/env nextflow

include { SCALE_DATA                            } from '../../../modules/local/scale_data'
include { RUN_PCA                               } from '../../../modules/local/run_pca'
include { QC_ELBOW_PLOT                         } from '../../../modules/local/qc_elbow_plot'
include { RUN_HARMONY                           } from '../../../modules/local/run_harmony'
include { RUN_UMAP                              } from '../../../modules/local/run_umap'
include { RUN_TSNE                              } from '../../../modules/local/run_tsne'
include { FIND_NEIGHBORS                        } from '../../../modules/local/find_neighbors'
include { FIND_CLUSTERS                         } from '../../../modules/local/find_clusters'
include { QC_DIM_PLOT_COUNTOUR as UMAP_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_DIM_PLOT_COUNTOUR as TSNE_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_HARMONY_PLOTS                      } from '../../../modules/local/qc_harmony_plots'
include { EXTRACT_SEURAT_CLUSTER_METADATA       } from '../../../modules/local/extract_seurat_cluster_metadata'
include { EXTRACT_SEURAT_REDUCED_DIMS           } from '../../../modules/local/extract_seurat_reduced_dims'
include { ADD_HARMONY_CLUSTER_TO_SEURAT         } from '../../../modules/local/add_harmony_cluster_to_seurat'

workflow INTEGRATE_HARMONY {
    take:
        ch_merged_xenium_obj    // channel: merged xenium objects
        dim_list                // list: list of dimensions to evaluate
        res_list                // list: list of resolutions to evaluate
        skip_tsne_plot          // boolean: whether to skip TSNE plot generation
        marker_gene_list        // file: marker gene list
        vf_nfeatures            // integer: number of variable features to select

    main:
        ch_versions = Channel.empty()

        // MODULE: Scale the merged xenium object
        SCALE_DATA ( ch_merged_xenium_obj )

        // MODULE: Run PCA for the xenium objects
        RUN_PCA ( SCALE_DATA.out.scaled_xenium_obj )

        // MODULE: Generate elbow plot
        QC_ELBOW_PLOT ( RUN_PCA.out.pca_xenium_obj )

        // MODULE: Run Harmony
        RUN_HARMONY ( RUN_PCA.out.pca_xenium_obj )

        // MODULE: Generate UMAPs for Harmony
        RUN_UMAP (
            RUN_HARMONY.out.integrated_xenium_obj
                .combine( Channel.from(dim_list) )
        )

        ch_umap_obj = RUN_UMAP.out.umap_xenium_obj

        if (!skip_tsne_plot) {
            // MODULE: Generate TSNE for Harmony
            RUN_TSNE (
                RUN_UMAP.out.umap_xenium_obj
            )   
            ch_umap_obj = RUN_TSNE.out.tsne_xenium_obj
        }

        // MODULE: Find Neighbors for Harmony Integrated object
        FIND_NEIGHBORS (
            ch_umap_obj
        )

        // MODULE: Find Clusters
        FIND_CLUSTERS (
            FIND_NEIGHBORS.out.find_neighbors_xenium_obj
                .combine( Channel.from(res_list) )
        )

        //
        // MODULE: Generate a dim plot with contours for UMAP
        //
        UMAP_DIM_PLOT (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

        ch_tsne_dim_plot = Channel.empty()
        if (!skip_tsne_plot) {
            //
            // MODULE: Generate a dim plot with contours for TSNE
            //
            TSNE_DIM_PLOT (
                FIND_CLUSTERS.out.find_clusters_xenium_obj
            )

            ch_tsne_dim_plot = TSNE_DIM_PLOT.out.countour_dim_plot
        }

        //
        // MODULE: Generate violin plots
        //

        if (marker_gene_list) {
            QC_HARMONY_PLOTS (
                FIND_CLUSTERS.out.find_clusters_xenium_obj
                    .combine( Channel.from(marker_gene_list) )
            )
        }

        //
        // MODULE: Extract cluster metadata
        //
        EXTRACT_SEURAT_CLUSTER_METADATA (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

        //
        // MODULE: Extract reduced dims
        //
        //TODO (see issue #22): check outputs for multiple resolutions (is the redudancy also happening here?)
        EXTRACT_SEURAT_REDUCED_DIMS (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

        //
        // MODULE: Add Harmony cluster info to Seurat object
        //
        ADD_HARMONY_CLUSTER_TO_SEURAT (
            RUN_PCA.out.pca_xenium_obj
                .join (
                    EXTRACT_SEURAT_CLUSTER_METADATA.out.cluster_metadata
                    .groupTuple()
                    .map{ meta, cm_file_list -> [meta, cm_file_list.flatten()]}
                )
                .join (
                    EXTRACT_SEURAT_REDUCED_DIMS.out.embeddings_csv
                        .groupTuple()
                        .map{ meta, e_file_list -> [meta, e_file_list.flatten()]}
                )
                .join (
                    EXTRACT_SEURAT_REDUCED_DIMS.out.loadings_csv
                        .groupTuple()
                        .map{ meta, l_file_list -> [meta, l_file_list.flatten()]}
                )
                .join (
                    EXTRACT_SEURAT_REDUCED_DIMS.out.stdev_csv
                        .groupTuple()
                        .map{ meta, s_file_list -> [meta, s_file_list.flatten()]}
                    )
        )

    emit:
        versions              = ch_versions
        integrated_xenium_obj = FIND_CLUSTERS.out.find_clusters_xenium_obj
        umap_dim_plot         = UMAP_DIM_PLOT.out.countour_dim_plot
        tsne_dim_plot         = ch_tsne_dim_plot

}
