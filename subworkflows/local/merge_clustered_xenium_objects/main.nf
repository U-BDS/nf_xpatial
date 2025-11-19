#!/usr/bin/env nextflow

include { EXTRACT_CLUSTER_METADATA   } from '../../../modules/local/extract_cluster_metadata'
include { EXTRACT_REDUCED_DIMS       } from '../../../modules/local/extract_reduced_dims'
include { MERGE_CSV                  } from '../../../modules/local/merge_csv'
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
        //TODO (see issue #22): check outputs for multiple resolutions (is the redudancy also happening here?)
        EXTRACT_REDUCED_DIMS (
            ch_clustered_xenium_obj
                .map{ meta, xenium_obj ->
                    def new_meta = [
                        id: meta.id,
                        normalization: meta.normalization,
                        dim: meta.dim ? meta.dim : meta.nPCs,
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
        
        ch_cluster_csvs.harmony.view()

        // MERGE_CSV (
        //     ch_cluster_csvs.banksy
        //         .groupTuple()
        // )

        //
        // MODULE: Connect BANKSY clusters
        //

        //
        // MODULE: Extract BANKSY connected clusters
        //

        //
        // MODULE: Add Harmony cluster info to Seurat object
        //
        // ADD_CLUSTER_DATA_TO_SEURAT (
        //     ch_merged_xenium_obj
        //         .join (
        //             EXTRACT_SEURAT_CLUSTER_METADATA.out.cluster_metadata
        //             .groupTuple()
        //             .map{ meta, cm_file_list -> 
        //                 def new_meta = [id: meta.id, normalization: meta.normalization]
        //                 [new_meta, cm_file_list.flatten()]}
        //         )
        //         .join (
        //             EXTRACT_SEURAT_REDUCED_DIMS.out.embeddings_csv
        //                 .groupTuple()
        //                 .map{ meta, e_file_list -> 
        //                     def new_meta = [id: meta.id, normalization: meta.normalization]
        //                     [new_meta, e_file_list.flatten()]}
        //         )
        //         .join (
        //             EXTRACT_SEURAT_REDUCED_DIMS.out.loadings_csv
        //                 .groupTuple()
        //                 .map{ meta, l_file_list -> 
        //                     def new_meta = [id: meta.id, normalization: meta.normalization]
        //                     [new_meta, l_file_list.flatten()]}
        //         )
        //         .join (
        //             EXTRACT_SEURAT_REDUCED_DIMS.out.stdev_csv
        //                 .groupTuple()
        //                 .map{ meta, s_file_list -> 
        //                     def new_meta = [id: meta.id, normalization: meta.normalization]
        //                     [new_meta, s_file_list.flatten()]}
        //             )
        // )

    emit:
        versions              = ch_versions

}
