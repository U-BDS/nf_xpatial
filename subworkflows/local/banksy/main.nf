#!/usr/bin/env nextflow

include { ADD_TISSUE_COORDS        } from '../../../modules/local/add_tissue_coords'
include { COMPILE_OBJECTS          } from '../../../modules/local/compile_objects'
include { MERGE_XENIUM_OBJECTS     } from '../../../modules/local/merge_xenium_objects'
include { CONVERT_SEURAT_TO_SPE    } from '../../../modules/local/convert_seurat_to_spe'
include { STAGGER_SPATIAL_COORDS   } from '../../../modules/local/stagger_spatial_coords'
include { CLUSTER_BANKSY           } from '../../../modules/local/cluster_banksy'
include { COMPUTE_BANKSY_MATRIX    } from '../../../modules/local/compute_banksy_matrix'
include { COMPUTE_BANKSY_PCA       } from '../../../modules/local/compute_banksy_pca'
include { RUN_HARMONY_BANKSY       } from '../../../modules/local/run_harmony_banksy'
include { RUN_UMAP_BANKSY          } from '../../../modules/local/run_umap_banksy'
include { EXTRACT_CLUSTER_METADATA } from '../../../modules/local/extract_cluster_metadata'
include { EXTRACT_PARAMS           } from '../../../modules/local/extract_params'
include { EXTRACT_XE_METADATA      } from '../../../modules/local/extract_xe_metadata'
include { EXTRACT_REDUCED_DIMS     } from '../../../modules/local/extract_reduced_dims'

workflow BANKSY {
    take:
        ch_comp_norm_xenium_obj // channel: compiled and normalized xenium objects
        lambda_list             // list: list of lambda values to evaluate
        k_geom_list             // list: list of k_geom values to evaluate
        nPCs_list               // list: list of nPCs values to evaluate
        res_list                // list: list of resolutions to evaluate

    main:
        ch_versions = Channel.empty()

        // MODULE: Add tissue coordiates to metadata
        ADD_TISSUE_COORDS ( ch_comp_norm_xenium_obj )

        // MODULE: Merge xenium objects
        MERGE_XENIUM_OBJECTS ( ADD_TISSUE_COORDS.out.tissue_coords_xenium_obj )

        // MODULE: Convert seurat object to spatial experiment object
        CONVERT_SEURAT_TO_SPE ( MERGE_XENIUM_OBJECTS.out.merged_xenium_obj )

        // MODULE: Stagger spatial coordinates
        STAGGER_SPATIAL_COORDS ( CONVERT_SEURAT_TO_SPE.out.spe_object )

        // MODULE: Compute BANKSY matrix
        COMPUTE_BANKSY_MATRIX (
            STAGGER_SPATIAL_COORDS.out.coord_staggered_spe_object
                .combine( Channel.from(k_geom_list) )
                .map {
                    meta, spe_obj, k_geom ->
                        meta.k_geom = k_geom
                        [meta, spe_obj, k_geom]
                }
        )

        // MODULE: Compute BANKSY PCAs
        COMPUTE_BANKSY_PCA (
            STAGGER_SPATIAL_COORDS.out.coord_staggered_spe_object
                .combine( Channel.from(lambda_list))
                .combine( Channel.from(nPCs_list) ).
                map {
                    meta, spe_obj, lambda, nPCs ->
                        meta.lambda = lambda
                        meta.nPCs = nPCs
                        [meta, spe_obj, lambda, nPCs]
                }
        )

        // MODULE: Run BANKSY Harmony
        RUN_HARMONY_BANKSY (
            COMPUTE_BANKSY_PCA.out.banksy_pca_spe_obj
                .map {
                    meta, spe_obj ->
                        [meta, spe_obj, meta.nPCs]
                }
        )

        // MODULE: Run BANKSY UMAP
        RUN_UMAP_BANKSY (
            RUN_HARMONY_BANKSY.out.banksy_pca_harmony_obj
        )

        // MODULE: CLUSTER BANKSY
        CLUSTER_BANKSY (
            RUN_UMAP_BANKSY.out.banksy_umap_spe_obj
                .combine( Channel.from(res_list) )
                .map {
                    meta, spe_obj, res ->
                        meta.res = res
                    [meta, spe_obj, meta.nPCs, meta.lambda, meta.res]
                }
        )

        // MODULE: Extract cluster metadata
        EXTRACT_CLUSTER_METADATA (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
        )

        // MODULE: Extract param data
        EXTRACT_PARAMS (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
        )

        // MODULE: Extract Xenium Explorer metadata
        EXTRACT_XE_METADATA (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
        )

        // MODULE: Extract Reduced Dims
        EXTRACT_REDUCED_DIMS (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
        )

    emit:
        versions = ch_versions

        ch_integrated_xenium_obj = Channel.empty()


}
