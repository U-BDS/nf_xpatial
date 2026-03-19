#!/usr/bin/env nextflow
include { COMPILE_ORDERED_OBJECTS } from '../../../subworkflows/local/compile_ordered_objects'

include { QC_IMAGE_DIM_PLOT                        } from '../../../modules/local/qc_image_dim_plot'
include { QC_VLN_PLOT                              } from '../../../modules/local/qc_vln_plot'
include { QC_FEATURE_SCATTER_PLOT                  } from '../../../modules/local/qc_feature_scatter_plot'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NCOUNT   } from '../../../modules/local/qc_image_feature_plot'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NFEATURE } from '../../../modules/local/qc_image_feature_plot'

workflow GENERAL_QC {
    take:
        ch_xenium_obj              // channel: annotated xenium data with associated marker genes
        skip_general_dim_plot      // bool: whether to skip the image dim plot
        skip_general_vln_plot      // bool: whether to skip the violin plot
        skip_general_scatter_plot  // bool: whether to skip the feature scatter plot
        skip_general_nfeature_plot // bool: whether to skip the nFeature image feature plot
        skip_general_ncount_plot   // bool: whether to skip the nCount image feature plot

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Compile ordered objects
        //
        COMPILE_ORDERED_OBJECTS (
            ch_xenium_obj
        )
        ch_compiled_obj = COMPILE_ORDERED_OBJECTS.out.compiled_obj

        // Create individual channels for each process to avoid race condition error on resume
        ch_dim_plot_in     = ch_compiled_obj
        ch_vln_plot_in     = ch_compiled_obj
        ch_scatter_plot_in = ch_compiled_obj
        ch_ifp_nfeature_in = ch_compiled_obj
        ch_ifp_ncount_in   = ch_compiled_obj

        ch_general_dim_plot = Channel.empty()
        if (!skip_general_dim_plot) {
            //
            // MODULE: Create an initial Image Dim Plot
            //
            QC_IMAGE_DIM_PLOT (
                ch_dim_plot_in
            )

            ch_general_dim_plot = QC_IMAGE_DIM_PLOT.out.image_dim_plot
            ch_versions = ch_versions.mix(QC_IMAGE_DIM_PLOT.out.versions)
        }

        ch_general_vln_plot = Channel.empty()
        if (!skip_general_vln_plot) {
            //
            // MODULE: Create a violin plot for nFeature and nCount
            //
            QC_VLN_PLOT (
                ch_vln_plot_in
            )

            ch_general_vln_plot = QC_VLN_PLOT.out.vln_plot
            ch_versions = ch_versions.mix(QC_VLN_PLOT.out.versions)
        }

        ch_general_scatter_plot = Channel.empty()
        if (!skip_general_scatter_plot) {
            //
            // MODULE: Create a feature scatter plot for nCount and nFeature
            //
            QC_FEATURE_SCATTER_PLOT (
                ch_scatter_plot_in
            )
            ch_versions = ch_versions.mix(QC_FEATURE_SCATTER_PLOT.out.versions)
            ch_general_scatter_plot = QC_FEATURE_SCATTER_PLOT.out.feature_scatter_plot
        }

        ch_general_nfeature_plot = Channel.empty()
        if (!skip_general_nfeature_plot) {
            //
            // MODULE: Create an image feature plot for nFeature
            //
            QC_IFP_NFEATURE(
                ch_ifp_nfeature_in
            )

            ch_general_nfeature_plot = QC_IFP_NFEATURE.out.image_feature_plot
            ch_versions = ch_versions.mix(QC_IFP_NFEATURE.out.versions)
        }

        ch_general_ncount_plot = Channel.empty()
        if (!skip_general_ncount_plot) {
            //
            // MODULE: Create an image feature plot for nCount
            //
            QC_IFP_NCOUNT(
                ch_ifp_ncount_in
            )

            ch_general_ncount_plot = QC_IFP_NCOUNT.out.image_feature_plot
            ch_versions = ch_versions.mix(QC_IFP_NCOUNT.out.versions)
        }

    emit:
        image_dim_plot              = ch_general_dim_plot
        vln_plot                    = ch_general_vln_plot
        feature_scatter_plot        = ch_general_scatter_plot
        image_feature_plot_ncount   = ch_general_ncount_plot
        image_feature_plot_nfeature = ch_general_nfeature_plot
        versions                    = ch_versions

}
