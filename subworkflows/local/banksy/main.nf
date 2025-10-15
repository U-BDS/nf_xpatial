#!/usr/bin/env nextflow

include { ADD_TISSUE_COORDS               } from '../../../modules/local/add_tissue_coords'
include { COMPILE_OBJECTS                 } from '../../../modules/local/compile_objects'
include { MERGE_XENIUM_OBJECTS            } from '../../../modules/local/merge_xenium_objects'
include { FIND_VARIABLE_FEATURES          } from '../../../modules/local/find_variable_features'
include { CONVERT_SEURAT_TO_SPE           } from '../../../modules/local/convert_seurat_to_spe'
include { SUBSET_VARIABLE_FEATURES        } from '../../../modules/local/subset_variable_features'
include { STAGGER_SPATIAL_COORDS          } from '../../../modules/local/stagger_spatial_coords'
include { CLUSTER_BANKSY                  } from '../../../modules/local/cluster_banksy'
include { COMPUTE_BANKSY_MATRIX           } from '../../../modules/local/compute_banksy_matrix'
include { COMPUTE_BANKSY_PCA              } from '../../../modules/local/compute_banksy_pca'
include { RUN_HARMONY_BANKSY              } from '../../../modules/local/run_harmony_banksy'
include { RUN_UMAP_BANKSY                 } from '../../../modules/local/run_umap_banksy'
include { EXTRACT_BANKSY_CLUSTER_METADATA } from '../../../modules/local/extract_banksy_cluster_metadata'
include { EXTRACT_PARAMS                  } from '../../../modules/local/extract_params'
include { EXTRACT_XE_METADATA             } from '../../../modules/local/extract_xe_metadata'
include { EXTRACT_BANKSY_REDUCED_DIMS     } from '../../../modules/local/extract_banksy_reduced_dims'
include { MERGE_CSV                       } from '../../../modules/local/merge_csv'
include { ADD_BANKSY_TO_SEURAT            } from '../../../modules/local/add_banksy_to_seurat'
include { QC_BANKSY_PLOTS                 } from '../../../modules/local/qc_banksy_plots'

workflow BANKSY {
    take:
        ch_comp_norm_xenium_obj // channel: compiled and normalized xenium objects
        lambda_list             // list: list of lambda values to evaluate
        k_geom_list             // list: list of k_geom values to evaluate
        nPCs_list               // list: list of nPCs values to evaluate
        res_list                // list: list of resolutions to evaluate
        skip_banksy_vf_filter   // boolean: whether to skip filtering to variable features
        vf_nfeatures            // integer: number of variable features to select


    main:
        ch_versions = Channel.empty()

        // MODULE: Add tissue coordiates to metadata
        ADD_TISSUE_COORDS ( ch_comp_norm_xenium_obj )

        // MODULE: Merge xenium objects
        MERGE_XENIUM_OBJECTS ( ADD_TISSUE_COORDS.out.tissue_coords_xenium_obj )

        ch_merged_xenium_obj = Channel.empty()
        if (!skip_banksy_vf_filter) {
            // MODULE: Find Variable Features
            FIND_VARIABLE_FEATURES ( 
                MERGE_XENIUM_OBJECTS.out.merged_xenium_obj,
                vf_nfeatures 
            )

            ch_merged_xenium_obj = FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj

        } else {
            ch_merged_xenium_obj = MERGE_XENIUM_OBJECTS.out.merged_xenium_obj
        }

        // MODULE: Convert seurat object to spatial experiment object
        CONVERT_SEURAT_TO_SPE ( ch_merged_xenium_obj )

        ch_spatial_exp = Channel.empty()
        if (!skip_banksy_vf_filter) {
            // MODULE: Subset to Variable Features
            SUBSET_VARIABLE_FEATURES ( CONVERT_SEURAT_TO_SPE.out.spe_object )
          
            ch_spatial_exp = SUBSET_VARIABLE_FEATURES.out.vf_subset_xenium_obj

        } else {
            ch_spatial_exp = CONVERT_SEURAT_TO_SPE.out.spe_object
        }

        // MODULE: Stagger spatial coordinates
        STAGGER_SPATIAL_COORDS ( ch_spatial_exp )

        // MODULE: Compute BANKSY matrix
        // This is to fix a race condition that occurs randomly on resumes
        ch_geom_list = Channel.of(k_geom_list).flatten()

        COMPUTE_BANKSY_MATRIX (
            STAGGER_SPATIAL_COORDS.out.coord_staggered_spe_object
                .combine( ch_geom_list )
        )

        // MODULE: Compute BANKSY PCAs

        // This is to fix a race condition that occurs randomly on resumes
        ch_lambda_npc = Channel.of(lambda_list)
            .flatten()
            .combine(
                Channel.of(nPCs_list).flatten()
            )

        COMPUTE_BANKSY_PCA (
            COMPUTE_BANKSY_MATRIX.out.banksy_mtx_spe_obj
                .combine( ch_lambda_npc )
        )

        // MODULE: Run BANKSY Harmony
        RUN_HARMONY_BANKSY (
            COMPUTE_BANKSY_PCA.out.banksy_pca_spe_obj
        )

        // MODULE: Run BANKSY UMAP
        RUN_UMAP_BANKSY (
            RUN_HARMONY_BANKSY.out.banksy_pca_harmony_obj
        )

        // MODULE: CLUSTER BANKSY
        CLUSTER_BANKSY (
            RUN_UMAP_BANKSY.out.banksy_umap_spe_obj
                .combine( Channel.of(res_list).flatten() )
        )

        // MODULE: Extract cluster metadata
        EXTRACT_BANKSY_CLUSTER_METADATA (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
        )

        // MODULE: Extract param data
        EXTRACT_PARAMS (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
                .map {
                    meta, csv, k_geom, lambda, nPCs, res ->
                        [meta,  csv]
                }
        )

        // MODULE: Extract Xenium Explorer metadata
        EXTRACT_XE_METADATA (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
                .map {
                    meta, csv, k_geom, lambda, nPCs, res ->
                        [meta,  csv]
                }
        )

        // MODULE: Extract Reduced Dims
        EXTRACT_BANKSY_REDUCED_DIMS (
            CLUSTER_BANKSY.out.banksy_cluster_spe_obj
                .map {
                    meta, csv, k_geom, lambda, nPCs, res ->
                        [meta,  csv]
                }
        )

        // MODULE: Merge cluster tsvs
        MERGE_CSV (
            EXTRACT_BANKSY_CLUSTER_METADATA.out.cluster_metadata
                .groupTuple()
        )

        // MODULE: Add BANKSY clusters to Xenium Object
        ADD_BANKSY_TO_SEURAT (
            MERGE_XENIUM_OBJECTS.out.merged_xenium_obj
                .join (
                    MERGE_CSV.out.merged_cluster_csv
                        .map {
                            meta, csv ->
                                keys_to_remove = ['lambda', 'k_geom', 'nPCs', 'res']
                                [meta.findAll { k, v -> !(k in keys_to_remove) }, csv]
                        }
                )
        )

        // MODULE: Generate QC plots for BANKSY clusters
        QC_BANKSY_PLOTS (
            ADD_BANKSY_TO_SEURAT.out.banksy_xenium_obj
                .combine (
                    EXTRACT_BANKSY_REDUCED_DIMS.out.banksy_umap_csv
                            .map {
                                meta, csv ->
                                    keys_to_remove = ['lambda', 'k_geom', 'nPCs', 'res']
                                    [meta.findAll { k, v -> !(k in keys_to_remove) }, csv]
                            }
                , by: 0)
        )

    emit:
        versions = ch_versions

        ch_integrated_xenium_obj = Channel.empty()


}
