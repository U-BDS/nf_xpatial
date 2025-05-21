#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// MODULE: Loaded from modules/local/
//
include { CREATE_XENIUM_OBJ                       } from '../modules/local/create_xenium_object'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { MANUAL_ANNOTATIONS_QC } from '../subworkflows/local/manual_annotations_qc/main'
include { MARKER_GENE_PAIRS_QC  } from '../subworkflows/local/marker_gene_pairs_qc/main'
include { CELL_SHAPE_QC         } from '../subworkflows/local/cell_shape_qc/main'
include { GENERAL_QC            } from '../subworkflows/local/general_qc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NF_XENIUM_ANALYSIS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    // Grab the header from the metadata sheet so we have values to group by in plots
    /*ch_samplesheet
        .map {
            meta, xenium_input, metadata, manual_annotation ->
                metadata_vals = metadata.text.readLines().first().split(',')
                metadata_vals.collect { metadata_val -> [meta, metadata_val] }
        }
        .flatMap{
            meta_metadata_val ->
                meta_metadata_val
        }
        .set { ch_metadata_vals }*/

    //
    // MODULE: Read in xenium matrix and add metadata
    //
    CREATE_XENIUM_OBJ (
        ch_samplesheet
            .map{
                meta, xenium_input, metadata, manual_annotation, gene_marker_list ->
                [meta, xenium_input, metadata]
            }
    )

    ch_xenium_obj = CREATE_XENIUM_OBJ.out.xenium_obj
    ch_versions = ch_versions.mix(CREATE_XENIUM_OBJ.out.versions)

    //
    // SUBWORKFLOW: Add manual annotations
    //
    MANUAL_ANNOTATIONS_QC (
        ch_samplesheet,
        ch_xenium_obj
    )

    //
    // SUBWORKFLOW: Perform qc for all marker gene pairings
    //
    MARKER_GENE_PAIRS_QC (
        ch_samplesheet.map { 
            meta, xenium_input, metadata, manual_annotation, marker_gene_list ->
            [meta, marker_gene_list]
        }
        .join(MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj, by: 0)
    )

    //
    // SUBWORKFLOW: Perform qc on cell shapes
    //
    CELL_SHAPE_QC (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // SUBWORKFLOW: Perform basic QC on the xenium objects
    //
    GENERAL_QC (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf_xenium_analysis_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
