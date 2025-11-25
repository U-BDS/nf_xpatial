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
include { ADD_TISSUE_COORDS                        } from '../modules/local/add_tissue_coords'
include { COMPILE_OBJECTS                          } from '../modules/local/compile_objects'
include { MERGE_XENIUM_OBJECTS                     } from '../modules/local/merge_xenium_objects'
include { FIND_VARIABLE_FEATURES                   } from '../modules/local/find_variable_features'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { MANUAL_ANNOTATIONS_QC               } from '../subworkflows/local/manual_annotations_qc/main'
include { SPATIAL_QC as SPATIAL_QC_PREFILTER  } from '../subworkflows/local/spatial_qc/main'
include { SPATIAL_QC as SPATIAL_QC_POSTFILTER } from '../subworkflows/local/spatial_qc/main'
include { NORMALIZE_DATA                      } from '../subworkflows/local/normalize_data/main'
include { CLUSTER_HARMONY                     } from '../subworkflows/local/cluster_harmony/main'
include { BANKSY                              } from '../subworkflows/local/banksy/main'
include { MERGE_CLUSTERED_XENIUM_OBJECTS      } from '../subworkflows/local/merge_clustered_xenium_objects/main'
include { CLUSTER_QC                          } from '../subworkflows/local/cluster_qc/main' 

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
    // SUBWORKFLOW: Generate spatial QC plots before filtering
    //
    SPATIAL_QC_PREFILTER (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj,
        params.marker_gene_list
    )

    //
    // MODULE: Filter the xenium data
    //
    FILTER_XENIUM_OBJ (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    //
    // SUBWORKFLOW: Generate spatial QC plots after filtering
    //
    SPATIAL_QC_POSTFILTER (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj,
        params.marker_gene_list
    )

    //
    // SUBWORKFLOW: Normalize xenium objects
    //
    NORMALIZE_DATA (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj,
        params.normalization_method ? params.normalization_method.split(',').collect { it.trim() } : []
    )

    //
    // MODULE: Add tissue coordiates to metadata
    //
    ADD_TISSUE_COORDS ( NORMALIZE_DATA.out.compiled_norm_objects )

    //
    // MODULE: Merge xenium objects
    //
    MERGE_XENIUM_OBJECTS ( ADD_TISSUE_COORDS.out.tissue_coords_xenium_obj )

    //
    // MODULE: Find Variable Features
    //
    FIND_VARIABLE_FEATURES ( 
        MERGE_XENIUM_OBJECTS.out.merged_xenium_obj,
        params.vf_nfeatures 
    )

    //
    // SUBWORKFLOW: Perform harmony integration on xenium objects
    //

    // Use the user-provided start, stop, and range values to generate a list of dimensions and resolutions
    def dim_list = params.selected_dim < 0
        ? (0..<( ((params.dim_stop + params.dim_step) - params.dim_start) / params.dim_step )).collect {params.dim_start + it * params.dim_step}
        : params.selected_dim

    def res_list = params.selected_res < 0
        ? (0..<( ((params.res_stop + params.res_step) - params.res_start) / params.res_step )).collect { params.res_start + it * params.res_step }
        : params.selected_res

    CLUSTER_HARMONY (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj,
        dim_list,
        res_list,
        params.skip_tsne_plot,
        params.marker_gene_list
    )

    ch_cluster_params = CLUSTER_HARMONY.out.harmony_clustered_xenium_obj
        .map {meta, xenium_obj ->
            def norm_method = meta.normalization
            [norm_method, meta]
        }
        .distinct()

    //
    // SUBWORKFLOW: Perform BANKSY clustering on xenium objects
    //
    BANKSY (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj,
        params.lambda_BANKSY.split(',').collect { it as Float },
        params.k_geom_BANKSY.split(',').collect { it as Integer },
        params.nPCs_BANKSY.split(',').collect { it as Integer },
        params.res_BANKSY.split(',').collect { it as Float },
        params.use_agf_BANKSY,
        params.skip_banksy_vf_filter
    )

    ch_cluster_params = ch_cluster_params
        .mix (
            BANKSY.out.banksy_clustered_xenium_obj
                .map {meta, xenium_obj ->
                    def norm_method = meta.normalization
                    [norm_method, meta]
                }
                .distinct()
        )

    //
    // SUBWORKFLOW: Merge Harmony and BANKSY clustered xenium objects
    //
    MERGE_CLUSTERED_XENIUM_OBJECTS (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj,
        BANKSY.out.banksy_clustered_xenium_obj
            .mix( CLUSTER_HARMONY.out.harmony_clustered_xenium_obj )
    )

    //
    // SUBWORKFLOW: Generate clustering QC plots
    //
    CLUSTER_QC (
        MERGE_CLUSTERED_XENIUM_OBJECTS.out.cluster_merged_obj
            .map { meta, xenium_obj ->
                def norm_method = meta.normalization
                [norm_method, meta, xenium_obj]
            }
            .combine(ch_cluster_params, by:0)
            .map { norm_method, merged_meta, xenium_obj, param_meta ->
                [ param_meta, xenium_obj]
            },
        params.marker_gene_list
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
