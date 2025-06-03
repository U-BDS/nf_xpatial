#!/usr/bin/env nextflow

include { QC_HISTOGRAM_PLOT                             } from '../../../modules/local/qc_histogram_plot'
include { COMPILE_OBJECTS                               } from '../../../modules/local/compile_objects'
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
        COMPILE_OBJECTS(
            ch_xenium_data
                .map{
                    meta, xenium_obj -> [xenium_obj]
                }
                .collect()
                .map{
                    [ [ 'id': 'compiled_filtered' ], it ]
                }
        )

        //
        // MODULE: Histogram of cell area
        //
        QC_HISTOGRAM_PLOT (
            COMPILE_OBJECTS.out.compiled_obj
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
        versions              = ch_versions


}