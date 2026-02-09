#!/usr/bin/env nextflow

include { RUN_BANKSY } from '../../../modules/local/run_banksy'
include { RUN_PCA    } from '../../../modules/local/run_pca'

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

        //
        // MODULE: Run BANKSY
        //

        RUN_BANKSY (
            ch_merged_xenium_obj
                .combine (Channel.of(lambda_list).flatten())
                .combine (Channel.of(k_geom_list).flatten())
                .map { meta, xenium_obj, lambda, k_geom ->
                    def new_meta = meta + [lambda: lambda, k_geom: k_geom]
                    [new_meta, xenium_obj]
                }
        )

        //
        // MODULE: Run PCAs
        //
        // nPcs is at script level, need to modify module file to allow nPCs to be passed on meta
        // We need to grab the maximum nPC and pass that in

        //
        // MODULE: Run Harmony
        //

        //
        // MODULE: Find Neighbors
        //

        //
        // MODULE: Find Clusters
        //


    emit:
        versions = ch_versions

}
