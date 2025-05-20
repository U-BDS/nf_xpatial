#!/usr/bin/env nextflow

include { ADD_MANUAL_ANNOTATIONS } from '../../../modules/local/add_manual_annotations'

workflow MANUAL_ANNOTATIONS_QC {
    take:
        ch_samplesheet   // channel: samplesheet read in from --input
        ch_xenium_obj    // channel: xenium object generated from samplesheet inputs
        ch_metadata_vals // channel: metadata values from the metadata sheet

    main:

        // Separate the samples that have manual annotations
        ch_sep_objects = ch_samplesheet
            .join(ch_xenium_obj)
            .map {
                meta, xenium_input, metadata, manual_annotation, xenium_rds ->
                    [meta, xenium_rds, manual_annotation]
            }
            .branch {
                meta, xenium_rds, manual_annotation ->
                    with_annotation: manual_annotation
                    no_annotation: true
            }

        //
        // MODULE: Add manual annotations where possible
        //
        ADD_MANUAL_ANNOTATIONS(
            ch_sep_objects.with_annotation
        )

        ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj
            .combine(ch_metadata_vals, by: 0)

        emit:
            annotated_xenium_obj = ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj
}

