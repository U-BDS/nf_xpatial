#!/usr/bin/env nextflow
include { CELL_SHAPE_QC         } from '../../../subworkflows/local/cell_shape_qc/main'
include { GENERAL_QC            } from '../../../subworkflows/local/general_qc/main'
include { CELL_AREA_QC          } from '../../../subworkflows/local/cell_area_qc/main'

workflow SPATIAL_QC {
    take:
        ch_xenium_obj                    // channel: xenium objects
        skip_cell_shape_qc               // bool: whether to skip the cell shape qc plots
        skip_cell_shape_prop_plot        // bool: whether to skip the cell shape proportion qc plots
        skip_cell_segm_prop_plot         // bool: whether to skip the cell segmentation proportion qc plots
        skip_general_qc                  // bool: whether to skip the general qc plots
        ch_skip_general_dim_plot         // bool: whether to skip the image dim plot
        ch_skip_general_vln_plot         // bool: whether to skip the violin plot
        ch_skip_general_scatter_plot     // bool: whether to skip the feature scatter plot
        ch_skip_general_nfeature_plot    // bool: whether to skip the nFeature image feature plot
        ch_skip_general_ncount_plot      // bool: whether to skip the nCount image feature plot
        skip_cell_area_qc                // bool: whether to skip the cell area qc plots
        skip_area_histogram_plot         // bool: whether to skip the area histogram plot
        skip_area_box_plot               // bool: whether to skip the area box plot
        skip_area_overlap_histogram_plot // bool: whether to skip the area overlapping histogram plot

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Generate qc plots for cell shapes
        //
        ch_cell_segm_prop_plot = Channel.empty()
        ch_cell_shape_prop_plot = Channel.empty()
        if (!skip_cell_shape_qc) {
            CELL_SHAPE_QC (
                ch_xenium_obj,
                skip_cell_shape_prop_plot,
                skip_cell_segm_prop_plot
            )
            ch_versions = ch_versions.mix(CELL_SHAPE_QC.out.versions)

            ch_cell_segm_prop_plot  = CELL_SHAPE_QC.out.cell_segm_plot
            ch_cell_shape_prop_plot = CELL_SHAPE_QC.out.cell_shape_plot
        }

        //
        // SUBWORKFLOW: Generate basic QC plots for xenium objects
        //
        ch_dim_plot = Channel.empty()
        ch_vln_plot = Channel.empty()
        ch_feature_scatter_plot = Channel.empty()
        ch_ncount_plot = Channel.empty()
        ch_nfeature_plot = Channel.empty()

        if (!skip_general_qc) {
            GENERAL_QC (
                ch_xenium_obj,
                ch_skip_general_dim_plot,
                ch_skip_general_vln_plot,
                ch_skip_general_scatter_plot,
                ch_skip_general_nfeature_plot,
                ch_skip_general_ncount_plot
            )
            ch_versions = ch_versions.mix(GENERAL_QC.out.versions)

            ch_dim_plot             = GENERAL_QC.out.image_dim_plot
            ch_vln_plot             = GENERAL_QC.out.vln_plot
            ch_feature_scatter_plot = GENERAL_QC.out.feature_scatter_plot
            ch_ncount_plot          = GENERAL_QC.out.image_feature_plot_ncount
            ch_nfeature_plot        = GENERAL_QC.out.image_feature_plot_nfeature
        }

        //
        // SUBWORKFLOW: Generate QC plots for cell area
        //
        ch_area_histogram_plot = Channel.empty()
        ch_area_box_plot = Channel.empty()
        ch_area_overlapping_histogram_plot = Channel.empty()

        if (!skip_cell_area_qc) {
            CELL_AREA_QC (
                ch_xenium_obj,
                skip_area_histogram_plot,
                skip_area_box_plot,
                skip_area_overlap_histogram_plot
            )
            ch_versions = ch_versions.mix(CELL_AREA_QC.out.versions)

            ch_area_histogram_plot             = CELL_AREA_QC.out.cell_area_histogram_plot
            ch_area_box_plot                   = CELL_AREA_QC.out.cell_area_box_plot
            ch_area_overlapping_histogram_plot = CELL_AREA_QC.out.cell_area_overlapping_histogram_plot
        }

    emit:
        versions                    = ch_versions

        // CELL_SHAPE_QC
        cell_segm_plot = ch_cell_segm_prop_plot
        cell_shape_plot = ch_cell_shape_prop_plot

        // GENERAL_QC
        dim_plot             = ch_dim_plot
        vln_plot             = ch_vln_plot
        feature_scatter_plot = ch_feature_scatter_plot
        ncount_plot          = ch_ncount_plot
        nfeature_plot        = ch_nfeature_plot

        // CELL_AREA_QC
        area_histogram_plot             = ch_area_histogram_plot
        area_box_plot                   = ch_area_box_plot
        area_overlapping_histogram_plot = ch_area_overlapping_histogram_plot



}
