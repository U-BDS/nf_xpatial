#!/usr/bin/env nextflow

include { QC_DIM_PLOT_COUNTOUR as UMAP_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_SPLIT_CLUSTER_PLOTS                } from '../../../modules/local/qc_split_cluster_plots'

workflow CLUSTER_QC {
    take:
        ch_clustered_xenium_obj    // channel: xenium objects
        marker_gene_list           // list: marker genes for violin plots

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Generate a dim plot with contours for UMAP
        //
        UMAP_DIM_PLOT (
            ch_clustered_xenium_obj
        )

        // //
        // // MODULE: Generate a compiled set of cluster plots split by sample
        // //
        // QC_SPLIT_CLUSTER_PLOTS (
        //     ch_clustered_xenium_obj
        // )

    emit:
        versions                    = ch_versions


}
