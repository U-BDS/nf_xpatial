#!/usr/bin/env nextflow

include { QC_DIM_PLOT_COUNTOUR as UMAP_DIM_PLOT } from '../../../modules/local/qc_dim_plot_countour'
include { QC_SPLIT_CLUSTER_PLOTS                } from '../../../modules/local/qc_split_cluster_plots'
include { QC_MARKER_VLN_PLOT                    } from '../../../modules/local/qc_marker_vln_plot'
include { QC_MARKER_DOT_PLOT                    } from '../../../modules/local/qc_marker_dot_plot'

workflow CLUSTER_QC {
    take:
        ch_clustered_xenium_obj    // channel: xenium objects
        marker_gene_list           // list: marker genes for violin plots
        skip_cluster_umap_plot     // bool: whether to skip the UMAP dim plot with contours
        skip_split_cluster_plot    // bool: whether to skip the split cluster plots
        skip_cluster_vln_plot      // bool: whether to skip the marker violin plots
        skip_cluster_dot_plot      // bool: whether to skip the marker dot plots

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Generate a dim plot with contours for UMAP
        //
        if (!skip_cluster_umap_plot) {
            UMAP_DIM_PLOT (
                ch_clustered_xenium_obj
            )
        }

        //
        // MODULE: Generate a compiled set of cluster plots split by sample
        //
        if (!skip_split_cluster_plot) {
            QC_SPLIT_CLUSTER_PLOTS (
                ch_clustered_xenium_obj
            )
        }

        if (marker_gene_list) {
            //
            // MODULE: Generate vln plots
            //
            if (!skip_cluster_vln_plot) {
                QC_MARKER_VLN_PLOT (
                    ch_clustered_xenium_obj,
                    marker_gene_list
                )
            }

            //
            // MODULE: Generate dot plots
            //
            if (!skip_cluster_dot_plot) {
                QC_MARKER_DOT_PLOT (
                    ch_clustered_xenium_obj,
                    marker_gene_list
                )
            }
        }

    emit:
        versions                    = ch_versions


}
