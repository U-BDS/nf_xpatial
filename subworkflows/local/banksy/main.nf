#!/usr/bin/env nextflow

include { CONVERT_SEURAT_TO_SPE           } from '../../../modules/local/convert_seurat_to_spe'
include { STAGGER_SPATIAL_COORDS          } from '../../../modules/local/stagger_spatial_coords'
include { COMPUTE_BANKSY_MATRIX           } from '../../../modules/local/compute_banksy_matrix'
include { COMPUTE_BANKSY_PCA              } from '../../../modules/local/compute_banksy_pca'
include { RUN_HARMONY_BANKSY              } from '../../../modules/local/run_harmony_banksy'
include { RUN_UMAP_BANKSY                 } from '../../../modules/local/run_umap_banksy'
include { CLUSTER_BANKSY                  } from '../../../modules/local/cluster_banksy'

workflow BANKSY {
    take:
        ch_merged_xenium_obj    // channel: merged xenium objects
        lambda_list             // list: list of lambda values to evaluate
        k_geom_list             // list: list of k_geom values to evaluate
        nPCs_list               // list: list of nPCs values to evaluate
        res_list                // list: list of resolutions to evaluate
        use_agf_BANKSY          // boolean: whether to use AGF in BANKSY


    main:
        ch_versions = Channel.empty()

        // MODULE: Convert seurat object to spatial experiment object
        CONVERT_SEURAT_TO_SPE (
            ch_merged_xenium_obj
                .map { meta, xenium_obj ->
                    def agf_val = use_agf_BANKSY == true ? 1 : 0
                    def new_meta = meta + [agf: agf_val]
                    [new_meta, xenium_obj]
                }
        )

        // MODULE: Stagger spatial coordinates
        STAGGER_SPATIAL_COORDS ( CONVERT_SEURAT_TO_SPE.out.spe_object )

        // MODULE: Compute BANKSY matrix
        // This is to fix a race condition that occurs randomly on resumes
        ch_geom_list = Channel.of(k_geom_list).flatten()

        COMPUTE_BANKSY_MATRIX (
            STAGGER_SPATIAL_COORDS.out.coord_staggered_spe_object
                .combine( ch_geom_list )
                .map {
                    meta, spe, k_geom ->
                        def new_meta = meta + [k_geom: k_geom]
                        [new_meta, spe]
                }
        )

        // MODULE: Compute BANKSY PCAs

        // This is to fix a race condition that occurs randomly on resumes
        ch_lambda_npc = Channel.of(lambda_list)
            .flatten()
            .combine(
                Channel.of(nPCs_list).flatten().max()
            )

        COMPUTE_BANKSY_PCA (
            COMPUTE_BANKSY_MATRIX.out.banksy_mtx_spe_obj
                .combine( ch_lambda_npc )
                .map {
                    meta, spe, lambda, max_nPC ->
                        def new_meta = meta + [lambda: lambda, max_nPC: max_nPC]
                        [new_meta, spe]
                }
        )

        // MODULE: Run BANKSY Harmony
        RUN_HARMONY_BANKSY (
            COMPUTE_BANKSY_PCA.out.banksy_pca_spe_obj
                .combine ( Channel.of(nPCs_list).flatten() )
                .map {
                    meta, spe, nPCs ->
                        def new_meta = meta + [nPCs: nPCs]
                        [new_meta, spe]
                }
        )

        // MODULE: Run BANKSY UMAP
        RUN_UMAP_BANKSY (
            RUN_HARMONY_BANKSY.out.banksy_pca_harmony_obj
        )

        // MODULE: CLUSTER BANKSY
        CLUSTER_BANKSY (
            RUN_UMAP_BANKSY.out.banksy_umap_spe_obj
                .combine( Channel.of(res_list).flatten() )
                .map {
                    meta, spe, res ->
                        def new_meta = meta + [res: res]
                        [new_meta, spe]
                }
        )

        ch_clustered_xenium_obj = CLUSTER_BANKSY.out.banksy_cluster_spe_obj
            .map { meta, xenium_obj ->
                def new_meta = meta + [clustering_method: 'BANKSY']
                [new_meta, xenium_obj]
            }

    emit:
        versions = ch_versions

        banksy_clustered_xenium_obj = ch_clustered_xenium_obj

}
