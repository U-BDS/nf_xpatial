#!/usr/bin/env nextflow

include { COMPILE_ORDERED_OBJECTS } from '../../../subworkflows/local/compile_ordered_objects'

include { QC_HISTOGRAM_PLOT                             } from '../../../modules/local/qc_histogram_plot'
include { COMPILE_DATAFRAMES as COMPILE_OUTLIER_DF      } from '../../../modules/local/compile_dataframes'
include { COMPILE_DATAFRAMES as COMPILE_OUTLIER_STAT_DF } from '../../../modules/local/compile_dataframes'
include { CALCULATE_OUTLIERS                            } from '../../../modules/local/calculate_outliers'
include { QC_BOX_PLOT                                   } from '../../../modules/local/qc_box_plot'
include { QC_OVERLAPPING_HISTOGRAM_PLOT                 } from '../../../modules/local/qc_overlapping_histogram_plot'

workflow CELL_AREA_QC {
    take:
        ch_xenium_data                   // channel: annotated xenium data
        skip_area_histogram_plot         // bool: whether to skip the area histogram plot
        skip_area_box_plot               // bool: whether to skip the area box plot
        skip_area_overlap_histogram_plot // bool: whether to skip the area overlapping histogram plot

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Compile objects
        //
        COMPILE_ORDERED_OBJECTS (
            ch_xenium_data
        )

        //
        // MODULE: Histogram of cell area
        //
        ch_area_histogram_plot = Channel.empty()
        if (!skip_area_histogram_plot) {
            QC_HISTOGRAM_PLOT (
                COMPILE_ORDERED_OBJECTS.out.compiled_obj
            )

            ch_area_histogram_plot = QC_HISTOGRAM_PLOT.out.histogram_plot
            ch_versions = ch_versions.mix(QC_HISTOGRAM_PLOT.out.versions)
        }

        //
        // MODULE: Calculate cell area outliers
        //
        CALCULATE_OUTLIERS (
            ch_xenium_data
        )
        ch_versions = ch_versions.mix(
            CALCULATE_OUTLIERS.out.versions
        )

        //
        // MODULE: Compile outlier csvs
        //
        COMPILE_OUTLIER_DF (
            CALCULATE_OUTLIERS.out.outlier_csv
                .map{
                    meta, outlier_csv -> [outlier_csv]
                }
                .collect()
                .map{
                    [ [ 'id': 'compiled_filtered' ], it ]
                }
        )

        ch_area_box_plot = Channel.empty()
        if (skip_area_box_plot) {
            //
            // MODULE: Box Plot of cell area with outliers
            //
            QC_BOX_PLOT (
                COMPILE_OUTLIER_DF.out.compiled_df
            )

            ch_area_box_plot = QC_BOX_PLOT.out.box_plot
        }

        ch_area_overlap_histogram_plot = Channel.empty()
        if (skip_area_overlap_histogram_plot) {
            //
            // MODULE: Overlapping Historgram Plot of cell area
            //
            QC_OVERLAPPING_HISTOGRAM_PLOT (
                COMPILE_OUTLIER_DF.out.compiled_df
            )

            ch_area_overlap_histogram_plot = QC_OVERLAPPING_HISTOGRAM_PLOT.out.overlapping_histogram_plot
        }


    emit:
        cell_area_obj                        = ch_xenium_data
        cell_area_histogram_plot             = ch_area_histogram_plot 
        cell_area_outliers                   = CALCULATE_OUTLIERS.out.outlier_csv
        cell_area_box_plot                   = ch_area_box_plot
        cell_area_overlapping_histogram_plot = ch_area_overlap_histogram_plot
        versions                             = ch_versions


}