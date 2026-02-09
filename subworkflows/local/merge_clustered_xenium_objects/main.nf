#!/usr/bin/env nextflow

include { EXTRACT_CLUSTER_METADATA   } from '../../../modules/local/extract_cluster_metadata'
include { EXTRACT_REDUCED_DIMS       } from '../../../modules/local/extract_reduced_dims'
include { MERGE_CSV                  } from '../../../modules/local/merge_csv'
include { CONNECT_CLUSTERS           } from '../../../modules/local/connect_clusters'
include { ADD_CLUSTER_DATA_TO_SEURAT } from '../../../modules/local/add_cluster_data_to_seurat'

workflow MERGE_CLUSTERED_XENIUM_OBJECTS {
    take:
        ch_merged_xenium_obj       // channel: pre-clustered merged xenium object
        ch_clustered_xenium_obj    // channel: post-clustered merged xenium objects

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Extract cluster metadata
        //
        EXTRACT_CLUSTER_METADATA (
            ch_clustered_xenium_obj
        )

        //
        // MODULE: Extract reduced dims
        //
        EXTRACT_REDUCED_DIMS (
            ch_clustered_xenium_obj
                .map{ meta, xenium_obj ->
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
                        dim: meta.dim ?: meta.nPCs,
                        clustering_method: meta.clustering_method
                    ]
                    [new_meta, xenium_obj]
                }
                .groupTuple()
                .map { meta, xenium_obj_list ->
                    [meta, xenium_obj_list.flatten().first()]
                }
        )
        //
        // MODULE: Merge cluster csvs
        //
        ch_cluster_csvs = EXTRACT_CLUSTER_METADATA.out.cluster_metadata
            .branch {
                meta, cluster_csv ->
                    banksy: meta.clustering_method == 'BANKSY'
                    harmony: meta.clustering_method == 'Harmony'
            }

        //
        // MODULE: Merge BANKSY cluster csvs
        //
        MERGE_CSV (
            ch_cluster_csvs.banksy
                .map{ meta, cluster_csv -> 
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
                    ]
                    [new_meta, cluster_csv]
                }
                .groupTuple()
        )

        //
        // MODULE: Connect BANKSY clusters
        //
        CONNECT_CLUSTERS (
            ch_merged_xenium_obj
                .join ( MERGE_CSV.out.merged_cluster_csv )
        )

        //
        // MODULE: Add Harmony cluster info to Seurat object
        //
        ADD_CLUSTER_DATA_TO_SEURAT (
            CONNECT_CLUSTERS.out.connected_xenium_obj
                .join (
                    ch_cluster_csvs.harmony
                    .map{ meta, cm_file_list -> 
                        def new_meta = [id: meta.id, normalization: meta.normalization]
                        [new_meta, cm_file_list]
                    }
                    .groupTuple()
                )
                .join (
                    EXTRACT_REDUCED_DIMS.out.embeddings_csv
                        .map { meta, e_file_list -> 
                            def new_meta = [id: meta.id, normalization: meta.normalization]
                            [new_meta, e_file_list]
                        }
                        .groupTuple()
                        .map{ meta, e_file_list -> 
                            [meta, e_file_list.flatten()]
                        }
                )
                .join (
                    EXTRACT_REDUCED_DIMS.out.loadings_csv
                        .map { meta, l_file_list -> 
                            def new_meta = [id: meta.id, normalization: meta.normalization]
                            [new_meta, l_file_list]
                        }
                        .groupTuple()
                        .map{ meta, l_file_list -> 
                            [meta, l_file_list.flatten()]
                        }
                )
                .join (
                    EXTRACT_REDUCED_DIMS.out.stdev_csv
                        .map { meta, s_file_list -> 
                            def new_meta = [id: meta.id, normalization: meta.normalization]
                            [new_meta, s_file_list]
                        }
                        .groupTuple()
                        .map{ meta, s_file_list -> 
                            [meta, s_file_list.flatten()]
                        }
                    )
        )

    emit:
        versions           = ch_versions
        cluster_merged_obj = ADD_CLUSTER_DATA_TO_SEURAT.out.all_cluster_xenium_obj

}
