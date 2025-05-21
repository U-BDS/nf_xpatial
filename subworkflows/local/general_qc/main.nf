#!/usr/bin/env nextflow
include { COMPILE_OBJECTS                          } from '../../../modules/local/compile_objects'
include { QC_IMAGE_DIM_PLOT                        } from '../../../modules/local/qc_image_dim_plot'
include { QC_VLN_PLOT                              } from '../../../modules/local/qc_vln_plot'
include { QC_FEATURE_SCATTER_PLOT                  } from '../../../modules/local/qc_feature_scatter_plot'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NCOUNT   } from '../../../modules/local/qc_image_feature_plot'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NFEATURE } from '../../../modules/local/qc_image_feature_plot'

workflow GENERAL_QC {
    take:
        ch_xenium_obj // channel: annotated xenium data with associated marker genes

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Compile objects into a list
        //
        COMPILE_OBJECTS (
            ch_xenium_obj
                .map{
                    meta, xenium_obj -> [xenium_obj]
                }
                .collect()
                .map{
                    [ [ 'id': 'compiled_RAW' ], it ]
                }
        )

        //
        // MODULE: Create an initial Image Dim Plot
        //
        QC_IMAGE_DIM_PLOT (
            COMPILE_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(QC_IMAGE_DIM_PLOT.out.versions)

        //
        // MODULE: Create a violin plot for nFeature and nCount
        //
        QC_VLN_PLOT (
            COMPILE_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(QC_VLN_PLOT.out.versions)

        //
        // MODULE: Create a feature scatter plot for nCount and nFeature
        //
        QC_FEATURE_SCATTER_PLOT (
            COMPILE_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(QC_FEATURE_SCATTER_PLOT.out.versions)

        //
        // MODULE: Create an image feature plot for nFeature
        //
        QC_IFP_NFEATURE(
            COMPILE_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(QC_IFP_NFEATURE.out.versions)

        //
        // MODULE: Create an image feature plot for nCount
        //
        QC_IFP_NCOUNT(
            COMPILE_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(QC_IFP_NCOUNT.out.versions)

    emit:
        image_dim_plot              = QC_IMAGE_DIM_PLOT.out.image_dim_plot
        vln_plot                    = QC_VLN_PLOT.out.vln_plot
        feature_scatter_plot        = QC_FEATURE_SCATTER_PLOT.out.feature_scatter_plot
        image_feature_plot_ncount   = QC_IFP_NCOUNT.out.image_feature_plot
        image_feature_plot_nfeature = QC_IFP_NFEATURE.out.image_feature_plot
        versions                    = ch_versions


}