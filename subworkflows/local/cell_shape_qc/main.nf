#!/usr/bin/env nextflow

include { COMPILE_ORDERED_OBJECTS } from '../../../subworkflows/local/compile_ordered_objects'

include { CLASSIFY_CELL_SHAPE                   } from '../../../modules/local/classify_cell_shape'
include { QC_PROPORTION_PLOT as CELL_SHAPE_PLOT } from '../../../modules/local/qc_proportion_plot'
include { QC_PROPORTION_PLOT as CELL_SEGM_PLOT  } from '../../../modules/local/qc_proportion_plot'

workflow CELL_SHAPE_QC {
    take:
        ch_xenium_data // channel: annotated xenium data with associated marker genes

    main:
        ch_versions = Channel.empty()
        //
        // MODULE: Classify the cell shape for all cells
        //
        CLASSIFY_CELL_SHAPE (
            ch_xenium_data
        )
        ch_versions = ch_versions.mix(CLASSIFY_CELL_SHAPE.out.versions)
    
        //
        // MODULE: Plot the cell shape classification 
        //
        // TODO: Need to obtain a user-provided list of cell ids, do NOT want to plot all cells
        // TODO: Is there a way to do this without adding a new column to the samplesheet?

        //
        // MODULE: Compile the cell shape dataframes into a single object
        //
        // TODO: Change this to use dataframes instead of objects
        COMPILE_ORDERED_OBJECTS (
            CLASSIFY_CELL_SHAPE.out.cell_shape_xenium_obj
        )
        //ch_versions = ch_versions.mix(COMPILE_OBJECTS.out.versions)

        //
        // MODULE: Plot the Cell Segmentation Distribution
        //
        CELL_SEGM_PLOT (
            COMPILE_ORDERED_OBJECTS.out.compiled_obj,
            "segmentation_method",
            "Sample"
        )
        ch_versions = ch_versions.mix(CELL_SEGM_PLOT.out.versions)

        //
        // MODULE: Plot the Cell Shape Distribution
        //
        CELL_SHAPE_PLOT (
            COMPILE_ORDERED_OBJECTS.out.compiled_obj,
            "shape_classification",
            "Sample"
        )
        ch_versions = ch_versions.mix(CELL_SHAPE_PLOT.out.versions)

    emit:
        cell_shape_xenium_obj = CLASSIFY_CELL_SHAPE.out.cell_shape_xenium_obj
        cell_shape_csv        = CLASSIFY_CELL_SHAPE.out.cell_shape_csv
        cell_shape_plot       = CELL_SHAPE_PLOT.out.proportion_plot
        cell_segm_plot        = CELL_SEGM_PLOT.out.proportion_plot
        versions              = ch_versions
}
