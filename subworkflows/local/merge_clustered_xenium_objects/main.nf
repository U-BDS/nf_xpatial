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

        ch_split_cluster_csvs = EXTRACT_CLUSTER_METADATA.out.cluster_metadata
            .branch {
                meta, cluster_csv ->
                    banksy: meta.clustering_method == 'BANKSY'
                    harmony: meta.clustering_method == 'Harmony'
                    banksy_seurat: meta.clustering_method == 'BANKSYSeurat'
            }

        //
        // MODULE: Extract reduced dims
        //

        EXTRACT_REDUCED_DIMS (
            ch_clustered_xenium_obj
                .map{ meta, xenium_obj ->
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
                        dim: meta.dim,
                        clustering_method: meta.clustering_method,
                        assay: meta.assay
                    ]
                    [new_meta, xenium_obj]
                }
                .groupTuple()
                .map { meta, xenium_obj_list ->
                    [meta, xenium_obj_list.flatten().first()]
                }
        )

        //
        // MODULE: Merge BANKSY cluster csvs
        //
        MERGE_CSV (
            ch_split_cluster_csvs.banksy
                .map{ meta, cluster_csv -> 
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
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
                .join ( MERGE_CSV.out.merged_cluster_csv )
        )

        // Create a channel of cluster csvs
        ch_split_cluster_csvs.harmony
            .mix( ch_split_cluster_csvs.banksy_seurat )
            .map { meta, cluster_csv -> 
                def new_meta = [
                    id: meta.id,
                    normalization: meta.normalization,
                    assay: meta.assay.replace('_BANKSY', '')
                ]
                [new_meta, cluster_csv]
            }
            .groupTuple()
            .set { ch_cluster_csvs }
        
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
                [meta, file_list.flatten()]
            }
            .multiMap { meta, files ->
                embeddings: files.findAll { it.name.contains('embeddings') } ? [meta, files.findAll { it.name.contains('embeddings') }] : null
                loadings: files.findAll { it.name.contains('loadings') } ? [meta, files.findAll { it.name.contains('loadings') }] : null  
                stdev: files.findAll { it.name.contains('stdev') } ? [meta, files.findAll { it.name.contains('stdev') }] : null
            }
            .set { ch_reduction_csvs }
        
        //
        // MODULE: Add Harmony cluster info to Seurat object
        //
        //CONNECT_CLUSTERS.out.connected_xenium_obj.view()

        // TODO: If the object already has cluster info, make sure to remove them
        ADD_CLUSTER_DATA_TO_SEURAT (
            CONNECT_CLUSTERS.out.connected_xenium_obj
                .join ( ch_cluster_csvs )
                .join ( ch_reduction_csvs.embeddings )
                .join ( ch_reduction_csvs.loadings )
                .join ( ch_reduction_csvs.stdev )
        )

    emit:
        versions           = ch_versions
        cluster_merged_obj = ADD_CLUSTER_DATA_TO_SEURAT.out.all_cluster_xenium_obj

}
