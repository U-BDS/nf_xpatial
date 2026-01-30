#!/usr/bin/env nextflow

include { SCALE_DATA                            } from '../../../modules/local/scale_data'
include { RUN_PCA                               } from '../../../modules/local/run_pca'
include { QC_ELBOW_PLOT                         } from '../../../modules/local/qc_elbow_plot'
include { RUN_HARMONY                           } from '../../../modules/local/run_harmony'
include { RUN_UMAP                              } from '../../../modules/local/run_umap'
include { RUN_TSNE                              } from '../../../modules/local/run_tsne'
include { FIND_NEIGHBORS                        } from '../../../modules/local/find_neighbors'
include { FIND_CLUSTERS                         } from '../../../modules/local/find_clusters'

workflow CLUSTER_HARMONY {
    take:
        ch_merged_xenium_obj    // channel: merged xenium objects
        dim_list                // list: list of dimensions to evaluate
        res_list                // list: list of resolutions to evaluate
        skip_tsne_plot          // boolean: whether to skip TSNE plot generation

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
                .map { meta, xenium_obj, dim ->
                    def new_meta = meta + [dim: dim]
                    return [new_meta, xenium_obj]
                }
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
                .map { meta, xenium_obj, res ->
                    def new_meta = meta + [res: res]
                    return [new_meta, xenium_obj]
                }
        )

        ch_clustered_xenium_obj = FIND_CLUSTERS.out.find_clusters_xenium_obj
            .map { meta, xenium_obj ->
                def new_meta = meta + [clustering_method: 'Harmony']
                [new_meta, xenium_obj]
            }

    emit:
        versions                     = ch_versions
        harmony_clustered_xenium_obj = ch_clustered_xenium_obj

}
