#!/usr/bin/env nextflow

include { ADD_MANUAL_ANNOTATIONS                              } from '../../../modules/local/add_manual_annotations'
include { QC_IMAGE_DIM_PLOT as MANUAL_ANNOTATION_IMG_DIM_PLOT } from '../../../modules/local/qc_image_dim_plot'
include { COMPILE_OBJECTS as COMPILE_MANUAL_ANNOTATIONS       } from '../../../modules/local/compile_objects'

workflow MANUAL_ANNOTATIONS_QC {
    take:
        ch_samplesheet   // channel: samplesheet read in from --input
        ch_xenium_obj    // channel: xenium object generated from samplesheet inputs

    main:
        ch_versions = Channel.empty()

        // Separate the samples that have manual annotations
        ch_sep_objects = ch_samplesheet
            .join(ch_xenium_obj)
            .map {
                meta, xenium_input, metadata, manual_annotation, gene_marker_list, xenium_rds ->
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
        ADD_MANUAL_ANNOTATIONS (
            ch_sep_objects.with_annotation
        )
        ch_versions = ch_versions.mix(ADD_MANUAL_ANNOTATIONS.out.versions)
        //ch_annotated_xenium_obj = ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj

        // Merge the xenium objects back together, removing the manual annotations from the channel
        ch_annotated_xenium_obj = ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj
        .mix (ch_sep_objects.no_annotation
            .map{
                meta, xenium_rds, manual_annotation ->
                    [meta, xenium_rds]
            }
        )

        //
        // MODULE: Merge the manual_annotations
        //
        COMPILE_MANUAL_ANNOTATIONS (
            ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj
                .map{
                    meta, xenium_obj -> [xenium_obj]
                }
                .collect()
                .map{
                    [ [ 'id': 'compiled_annotated' ], it ]
                }
        )
        //ch_versions = ch_versions.mix(COMPILE_MANUAL_ANNOTATIONS.out.versions)

        //
        // MODULE: Plot the manual annotations
        //
        MANUAL_ANNOTATION_IMG_DIM_PLOT (
            COMPILE_MANUAL_ANNOTATIONS.out.compiled_obj
        )
        ch_versions = ch_versions.mix(MANUAL_ANNOTATION_IMG_DIM_PLOT.out.versions)

        emit:
            annotated_xenium_obj     = ch_annotated_xenium_obj
            annotated_image_dim_plot = MANUAL_ANNOTATION_IMG_DIM_PLOT.out.image_dim_plot
            versions                 = ch_versions
}

