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
        ch_xenium_data // channel: annotated xenium data

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Calculate cell area
        //

        //
        // MODULE: Compile objects
        //
        COMPILE_ORDERED_OBJECTS (
            ch_xenium_data
        )

        //
        // MODULE: Histogram of cell area
        //
        QC_HISTOGRAM_PLOT (
            COMPILE_ORDERED_OBJECTS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(
            QC_HISTOGRAM_PLOT.out.versions
        )

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

        //
        // MODULE: Box Plot of cell area with outliers
        //
        QC_BOX_PLOT (
            COMPILE_OUTLIER_DF.out.compiled_df
        )

        //
        // MODULE: Overlapping Historgram Plot of cell area
        //
        QC_OVERLAPPING_HISTOGRAM_PLOT (
            COMPILE_OUTLIER_DF.out.compiled_df
        )


    emit:
        cell_area_obj                        = ch_xenium_data
        cell_area_outliers                   = CALCULATE_OUTLIERS.out.outlier_csv
        cell_area_box_plot                   = QC_BOX_PLOT.out.box_plot
        cell_area_overlapping_histogram_plot = QC_OVERLAPPING_HISTOGRAM_PLOT.out.overlapping_histogram_plot
        cell_area_histogram_plot             = QC_HISTOGRAM_PLOT.out.histogram_plot 
        versions                             = ch_versions


}