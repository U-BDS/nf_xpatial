#!/usr/bin/env nextflow

include { EXTRACT_CLUSTER_METADATA                                 } from '../../../modules/local/extract_cluster_metadata'
include { EXTRACT_CLUSTER_METADATA as EXTRACT_CONNECTED_CLUSTERS   } from '../../../modules/local/extract_cluster_metadata'
include { EXTRACT_REDUCED_DIMS                                     } from '../../../modules/local/extract_reduced_dims'
include { MERGE_CSV as MERGE_BANKSY_CSV                            } from '../../../modules/local/merge_csv'
include { MERGE_CSV as MERGE_CLUSTER_CSV                           } from '../../../modules/local/merge_csv'
include { CONNECT_CLUSTERS                                         } from '../../../modules/local/connect_clusters'
include { ADD_CLUSTER_DATA_TO_SEURAT                               } from '../../../modules/local/add_cluster_data_to_seurat'

workflow MERGE_CLUSTERED_XENIUM_OBJECTS {
    take:
        ch_merged_xenium_obj       // channel: pre-clustered merged xenium object
        ch_clustered_xenium_obj    // channel: post-clustered merged xenium objects

    main:
        ch_versions = Channel.empty()

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
                        clustering_method: meta.clustering_method,
                        assay: meta.assay
                    ]

                    if (meta.clustering_method == 'BANKSY' || meta.clustering_method == 'BANKSYSeurat') {
                        new_meta = new_meta + [
                            lambda: meta.lambda,
                            k_geom: meta.k_geom
                        ]
                    }

                    [new_meta, xenium_obj]
                }
                .groupTuple()
                .map { meta, xenium_obj_list ->
                    // Sort the list to ensure deterministic order
                    def sorted_list = xenium_obj_list.flatten().sort { it.toString() }
                    [meta, sorted_list[0]]
                }
        )

        //
        // MODULE: Extract cluster metadata
        //
        EXTRACT_CLUSTER_METADATA (
            ch_clustered_xenium_obj
        )

        ch_cluster_csvs = EXTRACT_CLUSTER_METADATA.out.cluster_metadata
            .branch {
                meta, cluster_csv ->
                    banksy: meta.clustering_method == 'BANKSY'
                    seurat: meta.clustering_method == 'Seurat'
                    banksy_seurat: meta.clustering_method == 'BANKSYSeurat'
            }


        //
        // MODULE: Merge BANKSY cluster csvs
        //
        MERGE_BANKSY_CSV (
            ch_cluster_csvs.banksy
                .map{ meta, cluster_csv -> 
                    def new_meta = [
                        id: meta.id,
                        clustering_method: meta.clustering_method,
                        normalization: meta.normalization,
                        lambda: meta.lambda,
                        assay: meta.assay
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
                .combine( 
                    MERGE_BANKSY_CSV.out.merged_cluster_csv
                        .map { meta, cluster_csv -> 
                            def new_meta = [
                                id: meta.id,
                                normalization: meta.normalization,
                                assay: meta.assay
                            ]
                            [new_meta, meta, cluster_csv]
                        }
                    , by: 0
                )
                .map { xenium_obj_meta, xenium_obj, cluster_csv_meta, cluster_csv ->
                    [cluster_csv_meta, xenium_obj, cluster_csv]
                }
        )

        //
        // MODULE: Extract connected cluster metadata
        //
        EXTRACT_CONNECTED_CLUSTERS (
            CONNECT_CLUSTERS.out.connected_xenium_obj
                .map { meta, xenium_obj -> 
                    def join_key = meta.normalization + meta.lambda
                    [join_key, meta, xenium_obj]}
                .combine (
                    ch_cluster_csvs.banksy
                        .map{ meta, cluster_csv -> 
                        def join_key = meta.normalization + meta.lambda
                        [join_key, meta]}
                    , by: 0
                )
                .map { norm_method, xenium_obj_meta, xenium_obj, cluster_csv_meta ->
                    [cluster_csv_meta, xenium_obj]
                }
        )

        //
        // MODULE: Merge all cluster csvs into a single csv
        //
        MERGE_CLUSTER_CSV (
            EXTRACT_CONNECTED_CLUSTERS.out.cluster_metadata
                .mix (ch_cluster_csvs.seurat)
                .mix (ch_cluster_csvs.banksy_seurat)
                .map { meta, cluster_csv -> 
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
                        assay: meta.assay.replace('_BANKSY', '')
                    ]
                    [new_meta, cluster_csv]
                }
                .groupTuple()
        )
        
        // Create and process the reduction csvs
        EXTRACT_REDUCED_DIMS.out.embeddings_csv
            .mix(EXTRACT_REDUCED_DIMS.out.loadings_csv)
            .mix(EXTRACT_REDUCED_DIMS.out.stdev_csv)
            .map { meta, file_list -> 
                def new_meta = [
                    id: meta.id,
                    normalization: meta.normalization,
                    assay: meta.assay.replace('_BANKSY', '')
                ]
                [new_meta, file_list]
            }
            .groupTuple()
            .map { meta, file_list -> 
                [meta, file_list.flatten().sort { it.toString() }]
            }
            .multiMap { meta, files ->
                embeddings: files.findAll { it.name.contains('embeddings') } ? [meta, files.findAll { it.name.contains('embeddings') }] : null
                loadings: files.findAll { it.name.contains('loadings') } ? [meta, files.findAll { it.name.contains('loadings') }] : null  
                stdev: files.findAll { it.name.contains('stdev') } ? [meta, files.findAll { it.name.contains('stdev') }] : null
            }
            .set { ch_reduction_csvs }
        
        //
        // MODULE: Add all cluster and dimension reductions to Seurat object
        //

        ADD_CLUSTER_DATA_TO_SEURAT (
            ch_merged_xenium_obj
                .join (
                    MERGE_CLUSTER_CSV.out.merged_cluster_csv
                )
                .join (
                    ch_reduction_csvs.embeddings
                )
                .join (
                    ch_reduction_csvs.loadings
                )
                .join (
                    ch_reduction_csvs.stdev
                )
        )

    emit:
        versions           = ch_versions
        cluster_merged_obj = ADD_CLUSTER_DATA_TO_SEURAT.out.all_cluster_xenium_obj

}
