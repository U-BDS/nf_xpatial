#!/usr/bin/env nextflow
include { MARKER_GENE_PAIRS_QC  } from '../../../subworkflows/local/marker_gene_pairs_qc/main'
include { CELL_SHAPE_QC         } from '../../../subworkflows/local/cell_shape_qc/main'
include { GENERAL_QC            } from '../../../subworkflows/local/general_qc/main'
include { CELL_AREA_QC          } from '../../../subworkflows/local/cell_area_qc/main'

workflow SPATIAL_QC {
    take:
        ch_xenium_obj    // channel: xenium objects
        marker_gene_list // file: marker gene list

    main:
        ch_versions = Channel.empty()

        //
        // SUBWORKFLOW: Generate qc plots for all marker gene pairings
        //
        MARKER_GENE_PAIRS_QC (
          ch_xenium_obj,
          marker_gene_list
        )
        ch_versions = ch_versions.mix(MARKER_GENE_PAIRS_QC.out.versions)

        //
        // SUBWORKFLOW: Generate qc plots for cell shapes
        //
        CELL_SHAPE_QC (
            ch_xenium_obj
        )
        ch_versions = ch_versions.mix(CELL_SHAPE_QC.out.versions)

        //
        // SUBWORKFLOW: Generate basic QC plots for xenium objects
        //
        GENERAL_QC (
            ch_xenium_obj
        )
        ch_versions = ch_versions.mix(GENERAL_QC.out.versions)

        //
        // SUBWORKFLOW: Generate QC plots for cell area
        //
        CELL_AREA_QC (
            ch_xenium_obj
        )
        ch_versions = ch_versions.mix(CELL_AREA_QC.out.versions)

    emit:
        versions                    = ch_versions


}
