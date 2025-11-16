#!/usr/bin/env nextflow

include { EXTRACT_SEURAT_CLUSTER_METADATA } from '../../../modules/local/extract_seurat_cluster_metadata'
include { EXTRACT_SEURAT_REDUCED_DIMS     } from '../../../modules/local/extract_seurat_reduced_dims'
include { ADD_CLUSTER_DATA_TO_SEURAT      } from '../../../modules/local/add_cluster_data_to_seurat'

workflow MERGE_CLUSTERED_XENIUM_OBJECTS {
    take:
        ch_merged_xenium_obj       // channel: pre-clustered merged xenium object
        ch_clustered_xenium_obj    // channel: post-clustered merged xenium objects

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Extract cluster metadata
        //
        EXTRACT_SEURAT_CLUSTER_METADATA (
            ch_clustered_xenium_obj
        )

        //
        // MODULE: Extract reduced dims
        //
        //TODO (see issue #22): check outputs for multiple resolutions (is the redudancy also happening here?)
        EXTRACT_SEURAT_REDUCED_DIMS (
            ch_clustered_xenium_obj
        )

            //
    // MODULE: Add Harmony cluster info to Seurat object
    //
    ADD_CLUSTER_DATA_TO_SEURAT (
        ch_merged_xenium_obj
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

}
