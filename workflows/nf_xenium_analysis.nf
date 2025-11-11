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
include { ADD_HARMONY_CLUSTER_TO_SEURAT            } from '../modules/local/add_harmony_cluster_to_seurat'
include { ADD_BANKSY_TO_SEURAT                     } from '../modules/local/add_banksy_to_seurat'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { MANUAL_ANNOTATIONS_QC               } from '../subworkflows/local/manual_annotations_qc/main'
include { SPATIAL_QC as SPATIAL_QC_PREFILTER  } from '../subworkflows/local/spatial_qc/main'
include { SPATIAL_QC as SPATIAL_QC_POSTFILTER } from '../subworkflows/local/spatial_qc/main'
include { NORMALIZE_DATA                      } from '../subworkflows/local/normalize_data/main'
include { INTEGRATE_HARMONY                   } from '../subworkflows/local/integrate_harmony/main'
include { BANKSY                              } from '../subworkflows/local/banksy/main'

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

    INTEGRATE_HARMONY (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj,
        dim_list,
        res_list,
        params.skip_tsne_plot,
        params.marker_gene_list
    )

    //
    // MODULE: Add Harmony cluster info to Seurat object
    //
    ADD_HARMONY_CLUSTER_TO_SEURAT (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj
            .join (
                INTEGRATE_HARMONY.out.harmony_cluster_metadata
                .groupTuple()
                .map{ meta, cm_file_list -> [meta, cm_file_list.flatten()]}
            )
            .join (
                INTEGRATE_HARMONY.out.harmony_embeddings
                    .groupTuple()
                    .map{ meta, e_file_list -> [meta, e_file_list.flatten()]}
            )
            .join (
                INTEGRATE_HARMONY.out.harmony_loadings
                    .groupTuple()
                    .map{ meta, l_file_list -> [meta, l_file_list.flatten()]}
            )
            .join (
                INTEGRATE_HARMONY.out.harmony_stdev
                    .groupTuple()
                    .map{ meta, s_file_list -> [meta, s_file_list.flatten()]}
                )
    )

    //
    // SUBWORKFLOW: Perform BANKSY clustering on xenium objects
    //
    BANKSY (
        FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj,
        params.lambda_BANKSY.split(',').collect { it as Float },
        params.k_geom_BANKSY.split(',').collect { it as Integer },
        params.nPCs_BANKSY.split(',').collect { it as Integer },
        params.res_BANKSY.split(',').collect { it as Float },
        params.skip_banksy_vf_filter
    )

    //
    // MODULE: Add BANKSY clusters to Xenium Object
    //
    ADD_BANKSY_TO_SEURAT (
        ADD_HARMONY_CLUSTER_TO_SEURAT.out.harmony_cluster_xenium_obj
            .join (
                BANKSY.out.merged_cluster_metadata
                    .map {
                        meta, csv ->
                            keys_to_remove = ['lambda', 'k_geom', 'nPCs', 'res']
                            [meta.findAll { k, v -> !(k in keys_to_remove) }, csv]
                    }
            )
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
