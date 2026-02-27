#!/usr/bin/env nextflow
include { MARKER_GENE_PAIRS_QC  } from '../../../subworkflows/local/marker_gene_pairs_qc/main'
include { CELL_SHAPE_QC         } from '../../../subworkflows/local/cell_shape_qc/main'
include { GENERAL_QC            } from '../../../subworkflows/local/general_qc/main'
include { CELL_AREA_QC          } from '../../../subworkflows/local/cell_area_qc/main'

workflow SPATIAL_QC {
    take:
        ch_xenium_obj                    // channel: xenium objects
        marker_gene_list                 // file: marker gene list
        skip_gene_list_filtering         // bool: whether to filter gene list used for marker gene pair qc
        skip_marker_gene_pair_qc         // bool: whether to skip the marker gene pair qc plots
        skip_marker_barnyard_plot        // bool: whether to skip the barnyard qc plots for marker gene pair qc,
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
        // SUBWORKFLOW: Generate qc plots for all marker gene pairings
        //
        if (!skip_marker_gene_pair_qc) {
            MARKER_GENE_PAIRS_QC (
                ch_xenium_obj,
                marker_gene_list,
                skip_gene_list_filtering,
                skip_marker_barnyard_plot
            )
            ch_versions = ch_versions.mix(MARKER_GENE_PAIRS_QC.out.versions)
        }

        //
        // SUBWORKFLOW: Generate qc plots for cell shapes
        //
        if (!skip_cell_shape_qc) {
            CELL_SHAPE_QC (
                ch_xenium_obj,
                skip_cell_shape_prop_plot,
                skip_cell_segm_prop_plot
            )
            ch_versions = ch_versions.mix(CELL_SHAPE_QC.out.versions)
        }

        //
        // SUBWORKFLOW: Generate basic QC plots for xenium objects
        //
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
        }

        //
        // SUBWORKFLOW: Generate QC plots for cell area
        //
        if (!skip_cell_area_qc) {
            CELL_AREA_QC (
                ch_xenium_obj,
                skip_area_histogram_plot,
                skip_area_box_plot,
                skip_area_overlap_histogram_plot
            )
            ch_versions = ch_versions.mix(CELL_AREA_QC.out.versions)
        }

    emit:
        versions                    = ch_versions


}
