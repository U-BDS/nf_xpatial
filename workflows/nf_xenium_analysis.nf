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
include { CREATE_XENIUM_OBJ                        } from '../modules/local/create_xenium_object'
include { FILTER_XENIUM_OBJ                        } from '../modules/local/filter_xenium_object'
include { COMPILE_OBJECTS as COMPILE_FILTERED_OBJS } from '../modules/local/compile_objects'
include { QC_VLN_PLOT as POST_FILTERING_VLN_PLOT   } from '../modules/local/qc_vln_plot'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { MANUAL_ANNOTATIONS_QC } from '../subworkflows/local/manual_annotations_qc/main'
include { MARKER_GENE_PAIRS_QC  } from '../subworkflows/local/marker_gene_pairs_qc/main'
include { CELL_SHAPE_QC         } from '../subworkflows/local/cell_shape_qc/main'
include { GENERAL_QC            } from '../subworkflows/local/general_qc/main'
include { CELL_AREA_QC          } from '../subworkflows/local/cell_area_qc/main'
include { NORMALIZE_DATA        } from '../subworkflows/local/normalize_data/main'
include { INTEGRATE_HARMONY     } from '../subworkflows/local/integrate_harmony/main'
include { BANKSY                } from '../subworkflows/local/banksy/main'

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
                meta, xenium_input, metadata, manual_annotation ->
                [meta, xenium_input, metadata]
            }
    )

    ch_xenium_obj = CREATE_XENIUM_OBJ.out.xenium_obj
    ch_versions = ch_versions.mix(CREATE_XENIUM_OBJ.out.versions)

    //
    // SUBWORKFLOW: Add manual annotations and produce qc plots
    //
    MANUAL_ANNOTATIONS_QC (
        ch_samplesheet,
        ch_xenium_obj
    )

    //
    // SUBWORKFLOW: Generate qc plots for all marker gene pairings
    //
    if (!params.marker_gene_list) {
        MARKER_GENE_PAIRS_QC (
            MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
                .combine(params.marker_gene_list)
        )
    }

    //
    // SUBWORKFLOW: Generate qc plots for cell shapes
    //
    CELL_SHAPE_QC (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // SUBWORKFLOW: Generate basic QC plots for xenium objects
    //
    GENERAL_QC (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // MODULE: Filter the xenium data
    //
    FILTER_XENIUM_OBJ (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // SUBWORKFLOW: Generate QC plots for cell area
    //
    //TODO: Does this need to be run post-filtering?
    CELL_AREA_QC (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj
    )

    //
    // MODULE: Post filtering Violin Plots
    //
    COMPILE_FILTERED_OBJS (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj
            .map{
                meta, xenium_obj -> [xenium_obj]
            }
            .collect()
            .map{
                [ [ 'id': 'compiled_FILTERED' ], it ]
            }
    )

    //
    // MODULE: Post filtering Violin Plots
    //

    POST_FILTERING_VLN_PLOT (
        COMPILE_FILTERED_OBJS.out.compiled_obj
    )

    //
    // SUBWORKFLOW: Normalize xenium objects
    //
    NORMALIZE_DATA (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj,
        params.normalization_method ? params.normalization_method.split(',').collect { it.trim() } : []
    )

    //
    // SUBWORKFLOW: Perform harmony integration on xenium objects
    //

    // Use the user-provided start, stop, and range values to generate a list of dimensions and resolutions
    def dim_list = params.selected_dim < 0
        ? (0..<( (params.dim_stop - params.dim_start) / params.dim_step )).collect {params.dim_start + it * params.dim_step}
        : params.selected_dim

    def res_list = params.selected_res < 0
        ? (0..<( (params.res_stop - params.res_start) / params.res_step )).collect { params.res_start + it * params.res_step }
        : params.selected_res

    INTEGRATE_HARMONY (
        NORMALIZE_DATA.out.compiled_norm_objects,
        dim_list,
        res_list
    )

    //
    // SUBWORKFLOW: Perform BANKSY clustering on xenium objects
    //

    // Generate the list of banksy values from user-provided lists
    BANKSY (
        NORMALIZE_DATA.out.compiled_norm_objects,
        params.lambda_BANKSY.split(',').collect { it as Float },
        params.k_geom_BANKSY.split(',').collect { it as Integer },
        params.nPCs_BANKSY.split(',').collect { it as Integer },
        params.res_BANKSY.split(',').collect { it as Float },
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
