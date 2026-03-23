#!/usr/bin/env nextflow


include { SUBSET_VARIABLE_FEATURES        } from '../../../modules/local/subset_variable_features'
include { RUN_BANKSY     } from '../../../modules/local/run_banksy'
include { RUN_PCA        } from '../../../modules/local/run_pca'
include { RUN_UMAP       } from '../../../modules/local/run_umap'
include { RUN_HARMONY    } from '../../../modules/local/run_harmony'
include { FIND_NEIGHBORS } from '../../../modules/local/find_neighbors'
include { FIND_CLUSTERS  } from '../../../modules/local/find_clusters'


workflow CLUSTER_BANKSY_SEURAT_WRAPPER {
    take:
        ch_merged_xenium_obj    // channel: merged xenium objects
        lambda_list             // list: list of lambda values to evaluate
        k_geom_list             // list: list of k_geom values to evaluate
        nPCs_list               // list: list of nPCs values to evaluate
        res_list                // list: list of resolutions to evaluate
        use_agf_BANKSY          // boolean: whether to use AGF in BANKSY


    main:
        ch_versions = Channel.empty()

        // MODULE: Subset to Variable Features
        SUBSET_VARIABLE_FEATURES (
            ch_merged_xenium_obj
        )

        // MODULE: Run BANKSY
        RUN_BANKSY (
            SUBSET_VARIABLE_FEATURES.out.vf_subset_xenium_obj
                .combine (Channel.of(lambda_list).flatten())
                .combine (Channel.of(k_geom_list).flatten())
                .map { 
                    meta, xenium_obj, lambda, k_geom ->
                        def new_meta = meta + [lambda: lambda, k_geom: k_geom]
                        [new_meta, xenium_obj]
                }
        )

        ch_banksy_xenium_obj = RUN_BANKSY.out.banksy_xenium_obj
            .map {
                meta, xenium_obj ->
                    def new_meta = meta + [assay: "${meta.assay}_BANKSY"]
                    [new_meta, xenium_obj]
            }

        // MODULE: Run PCAs
        // TODO: nPcs is at script level, need to modify module file to allow nPCs to be passed on meta
        //      We need to grab the maximum nPC and pass that in
        // TODO: Do we need variable features in here? We run variable features before
        RUN_PCA (
            ch_banksy_xenium_obj
        )

        // // MODULE: Run Harmony
        RUN_HARMONY ( RUN_PCA.out.pca_xenium_obj )

        // MODULE: Generate UMAPs for Harmony
        RUN_UMAP (
            RUN_HARMONY.out.integrated_xenium_obj
                .combine( Channel.from(nPCs_list) )
                .map { meta, xenium_obj, dim ->
                    def new_meta = meta + [dim: dim]
                    return [new_meta, xenium_obj]
                }
        )

        // MODULE: Find Neighbors for Harmony Integrated object
        FIND_NEIGHBORS (
            RUN_UMAP.out.umap_xenium_obj
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
                def new_meta = meta + [clustering_method: 'BANKSYSeurat']
                [new_meta, xenium_obj]
            }


    emit:
        versions = ch_versions
        banksy_xenium_obj = ch_banksy_xenium_obj
        clustered_xenium_obj = ch_clustered_xenium_obj


}
