#!/usr/bin/env nextflow

include { MERGE_XENIUM_OBJECTS                 } from '../../../modules/local/merge_xenium_objects'
include { SCALE_DATA                           } from '../../../modules/local/scale_data'
include { RUN_PCA                              } from '../../../modules/local/run_pca'
include { QC_ELBOW_PLOT                        } from '../../../modules/local/qc_elbow_plot'
include { RUN_HARMONY                          } from '../../../modules/local/run_harmony'
include { RUN_UMAP                             } from '../../../modules/local/run_umap'
include { RUN_TSNE                             } from '../../../modules/local/run_tsne'
include { FIND_NEIGHBORS                       } from '../../../modules/local/find_neighbors'
include { FIND_CLUSTERS                        } from '../../../modules/local/find_clusters'
include { QC_DIM_PLOT_COUNTOUR as UMAP_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_DIM_PLOT_COUNTOUR as TSNE_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_VLN_PLOT                          } from '../../../modules/local/qc_vln_plot'
include { QC_IMAGE_DIM_PLOT                    } from '../../../modules/local/qc_image_dim_plot'

workflow INTEGRATE_HARMONY {
    take:
        ch_comp_norm_xenium_obj // channel: compiled and normalized xenium objects
        dim_list                // list: list of dimensions to evaluate
        res_list                // list: list of resolutions to evaluate

    main:
        ch_versions = Channel.empty()

        // MODULE: Merge xenium objects
        MERGE_XENIUM_OBJECTS ( ch_comp_norm_xenium_obj )

        // MODULE: Scale the merged xenium object
        SCALE_DATA ( MERGE_XENIUM_OBJECTS.out.merged_xenium_obj )

        // MODULE: Run PCA for the xenium objects
        RUN_PCA ( SCALE_DATA.out.scaled_xenium_obj )

        // MODULE: Generate elbow plot
        QC_ELBOW_PLOT ( RUN_PCA.out.pca_xenium_obj )

        // MODULE: Run Harmony
        RUN_HARMONY ( RUN_PCA.out.pca_xenium_obj )

        // MODULE: Generate UMAPs for Harmony
        RUN_UMAP (
            RUN_HARMONY.out.integrated_xenium_obj.combine( Channel.from(dim_list) )
        )

        // MODULE: Generate TSNE for Harmony
        RUN_TSNE (
            RUN_UMAP.out.umap_xenium_obj
        )

        // MODULE: Find Neighbors for Harmony Integrated object
        FIND_NEIGHBORS (
            RUN_TSNE.out.tsne_xenium_obj
        )

        // MODULE: Find Clusters
        FIND_CLUSTERS (
            FIND_NEIGHBORS.out.find_neighbors_xenium_obj
                .combine( Channel.from(res_list) )
                .map {
                    meta, xenium_obj, dim, res ->
                        meta.dim = dim
                        meta.res = res
                        [meta, xenium_obj, res]
                }
        )

        //
        // MODULE: Generate a dim plot with contours for UMAP
        //
        UMAP_DIM_PLOT (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

        //
        // MODULE: Generate a dim plot with contours for TSNE
        //
        TSNE_DIM_PLOT (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

        //
        // MODULE: Generate violin plots
        //
        QC_VLN_PLOT (
            FIND_CLUSTERS.out.find_clusters_xenium_obj
        )

    emit:
        versions = ch_versions

        ch_integrated_xenium_obj = Channel.empty()


}
