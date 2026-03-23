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
include { ADD_METADATA                             } from '../modules/local/add_metadata'
include { FILTER_XENIUM_OBJ                        } from '../modules/local/filter_xenium_object'
include { COMPILE_OBJECTS as COMPILE_FILTERED_OBJS } from '../modules/local/compile_objects'
include { ADD_TISSUE_COORDS                        } from '../modules/local/add_tissue_coords'
include { COMPILE_OBJECTS                          } from '../modules/local/compile_objects'
include { MERGE_XENIUM_OBJECTS                     } from '../modules/local/merge_xenium_objects'
include { COMPILE_IMAGES_TO_VIDEO                  } from '../modules/local/compile_images_to_video'
include { RENDER_SUMMARY_REPORT                    } from '../modules/local/render_summary_report'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//
include { MANUAL_ANNOTATIONS_QC               } from '../subworkflows/local/manual_annotations_qc/main'
include { SPATIAL_QC as SPATIAL_QC_PREFILTER  } from '../subworkflows/local/spatial_qc/main'
include { SPATIAL_QC as SPATIAL_QC_POSTFILTER } from '../subworkflows/local/spatial_qc/main'
include { NORMALIZE_DATA                      } from '../subworkflows/local/normalize_data/main'
include { GET_VARIABLE_FEATURES               } from '../subworkflows/local/get_variable_features'
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

    ch_input_types = ch_samplesheet
        .branch {
            meta, xenium_input, metadata, manual_annotation ->
                seurat_obj: xenium_input.endsWith('.rds')
                xenium_out: true
        }

    //
    // MODULE: Read in xenium matrix and add metadata
    //
    CREATE_XENIUM_OBJ (
        ch_input_types.xenium_out
            .map{
                meta, xenium_input, metadata, manual_annotation ->
                [meta, xenium_input]
            }
    )

    ch_xenium_obj = CREATE_XENIUM_OBJ.out.xenium_obj
    ch_versions = ch_versions.mix(CREATE_XENIUM_OBJ.out.versions)

    ADD_METADATA (
        ch_xenium_obj
            .join ( ch_samplesheet )
            .mix ( ch_input_types.seurat_obj )
            .map { meta, xenium_obj, xenium_input, metadata, manual_annotation ->
                [meta, xenium_obj, metadata]
            }
    )

    //
    // SUBWORKFLOW: Add manual annotations and produce qc plots
    //
    MANUAL_ANNOTATIONS_QC (
        ch_samplesheet,
        ADD_METADATA.out.metadata_xenium_obj,
        params.skip_man_ann_dim_plot
    )

    if (!params.skip_qc) {
        //
        // SUBWORKFLOW: Generate spatial QC plots before filtering
        //
        SPATIAL_QC_PREFILTER (
            MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj,
            params.marker_gene_list,
            params.skip_gene_list_filtering,
            params.skip_pre_marker_gene_qc,
            params.skip_pre_marker_barnyard_plot,
            params.skip_pre_cell_shape_qc,
            params.skip_pre_cell_shape_prop_plot,
            params.skip_pre_cell_segm_prop_plot,
            params.skip_pre_general_qc,
            params.skip_pre_general_dim_plot,
            params.skip_pre_general_vln_plot,
            params.skip_pre_general_scatter_plot,
            params.skip_pre_general_nfeature_plot,
            params.skip_pre_general_ncount_plot,
            params.skip_pre_cell_area_qc,
            params.skip_pre_area_histogram_plot,
            params.skip_pre_area_box_plot,
            params.skip_pre_area_overlap_histogram_plot
        )
    }

    //
    // MODULE: Filter the xenium data
    //
    FILTER_XENIUM_OBJ (
        MANUAL_ANNOTATIONS_QC.out.annotated_xenium_obj
    )

    if (!params.skip_qc) {
        //
        // SUBWORKFLOW: Generate spatial QC plots after filtering
        //
        SPATIAL_QC_POSTFILTER (
            FILTER_XENIUM_OBJ.out.filtered_xenium_obj,
            params.marker_gene_list,
            params.skip_gene_list_filtering,
            params.skip_post_marker_gene_qc,
            params.skip_post_marker_barnyard_plot,
            params.skip_post_cell_shape_qc,
            params.skip_post_cell_shape_prop_plot,
            params.skip_post_cell_segm_prop_plot,
            params.skip_post_general_qc,
            params.skip_post_general_dim_plot,
            params.skip_post_general_vln_plot,
            params.skip_post_general_scatter_plot,
            params.skip_post_general_nfeature_plot,
            params.skip_post_general_ncount_plot,
            params.skip_post_cell_area_qc,
            params.skip_post_area_histogram_plot,
            params.skip_post_area_box_plot,
            params.skip_post_area_overlap_histogram_plot
        )
    }

    //
    // SUBWORKFLOW: Normalize xenium objects
    //
    NORMALIZE_DATA (
        FILTER_XENIUM_OBJ.out.filtered_xenium_obj,
        params.normalization_method ? params.normalization_method.split(',').collect { it.trim() } : [],
        params.skip_qc || params.skip_norm_ncount,
        params.skip_qc || params.skip_norm_nfeature
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
    // SUBWORKFLOW: Find Variable Features
    //
    GET_VARIABLE_FEATURES ( 
        MERGE_XENIUM_OBJECTS.out.merged_xenium_obj,
        params.skip_banksy_vf_filter
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
        GET_VARIABLE_FEATURES.out.vf_xenium_obj,
        dim_list,
        res_list,
        params.skip_qc || params.skip_tsne_plot
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
        params.skip_banksy_vf_filter ? GET_VARIABLE_FEATURES.out.vf_xenium_obj : GET_VARIABLE_FEATURES.out.vf_subset_obj,
        params.lambda_BANKSY.split(',').collect { it as Float },
        params.k_geom_BANKSY.split(',').collect { it as Integer },
        params.nPCs_BANKSY.split(',').collect { it as Integer },
        params.res_BANKSY.split(',').collect { it as Float },
        params.use_agf_BANKSY
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
    // SUBWORKFLOW: Perform BANKSY clustering on xenium objects using Seurat wrapper
    //

    CLUSTER_BANKSY_SEURAT_WRAPPER (
        params.skip_banksy_vf_filter ? GET_VARIABLE_FEATURES.out.vf_xenium_obj : GET_VARIABLE_FEATURES.out.vf_subset_obj,
        params.lambda_BANKSY.split(',').collect { it as Float },
        params.k_geom_BANKSY.split(',').collect { it as Integer },
        params.nPCs_BANKSY.split(',').collect { it as Integer },
        params.res_BANKSY.split(',').collect { it as Float },
        params.use_agf_BANKSY
    )

    ch_cluster_params = CLUSTER_BANKSY_SEURAT_WRAPPER.out.clustered_xenium_obj
        .map {meta, xenium_obj ->
            def norm_method = meta.normalization
            [norm_method, meta]
        }
        .distinct()

    //
    // SUBWORKFLOW: Merge Harmony and BANKSY clustered xenium objects
    //
    MERGE_CLUSTERED_XENIUM_OBJECTS (
        GET_VARIABLE_FEATURES.out.vf_xenium_obj,
        BANKSY.out.banksy_clustered_xenium_obj
            .mix( CLUSTER_HARMONY.out.harmony_clustered_xenium_obj )
    )

    //
    // SUBWORKFLOW: Generate clustering QC plots
    //
    if (!params.skip_qc) {
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
            params.marker_gene_list ?: GET_VARIABLE_FEATURES.out.gene_list,
            params.skip_qc || params.skip_cluster_umap_plot,
            params.skip_qc || params.skip_cluster_split_plot,
            params.skip_qc || params.skip_cluster_vln_plot,
            params.skip_qc || params.skip_cluster_dot_plot
        )
    }

    //
    // MODULE: Compile images into video
    //
    COMPILE_IMAGES_TO_VIDEO (
        CLUSTER_QC.out.umap_plot
            .map { meta, umap_plot ->
                def new_meta = [ 
                    'id': meta.id,
                    'normalization': meta.normalization,
                    'clustering_method': meta.clustering_method,
                    'plot_type': 'umap'
                ]
                [new_meta, umap_plot]
            }
            .mix (
                CLUSTER_QC.out.split_cluster_plot
                    .map { meta, split_cluster_plot ->
                        def new_meta = [ 
                            'id': meta.id,
                            'normalization': meta.normalization,
                            'clustering_method': meta.clustering_method,
                            'plot_type': 'split_cluster'
                        ]
                        [new_meta, split_cluster_plot]
                    }
            )
            .mix (
                CLUSTER_QC.out.marker_vln_plot
                    .map { meta, marker_vln_plot ->
                        def new_meta = [ 
                            'id': meta.id,
                            'normalization': meta.normalization,
                            'clustering_method': meta.clustering_method,
                            'plot_type': 'vln_plot'
                        ]
                        [new_meta, marker_vln_plot]
                    }
            )
            .mix (
                CLUSTER_QC.out.marker_dot_plot
                    .map { meta, marker_dot_plot ->
                        def new_meta = [ 
                            'id': meta.id,
                            'normalization': meta.normalization,
                            'clustering_method': meta.clustering_method,
                            'plot_type': 'dot_plot'
                        ]
                        [new_meta, marker_dot_plot]
                    }
            )
            .groupTuple()
            .map { meta, plot_list -> [meta, plot_list.sort{it.toString()}] }
    )

    //
    // MODULE: Create the summary report for the analysis
    //
    RENDER_SUMMARY_REPORT(
    Channel.of(file("$baseDir/assets/report_template.Rmd"))
        .map { report_template -> 
            def meta = ['id': 'summary_report']
            [meta, report_template]
        }
        .combine ( Channel.of(file("$baseDir/assets/style.css")) )
        .combine ( Channel.of("${params.min_nCount}") )
        .combine ( Channel.of("${params.min_nFeature}") )
        .combine ( Channel.of("${dim_list}") )
        .combine ( Channel.of("${res_list}") )
        .combine ( Channel.of("${params.lambda_BANKSY}") )
        .combine ( Channel.of("${params.k_geom_BANKSY}") )
        .combine ( Channel.of("${params.nPCs_BANKSY}") )
        .combine ( Channel.of("${params.res_BANKSY}") )
        .combine (
            FILTER_XENIUM_OBJ.out.filtered_stats_csv
                .map{ meta, filtered_stat -> [1, filtered_stat] }
                .groupTuple()
                .map{ key, filtered_file_list -> [filtered_file_list] }
        )
        .combine ( SPATIAL_QC_POSTFILTER.out.dim_plot.map {meta, dim_plot -> [dim_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.vln_plot.map {meta, vln_plot -> [vln_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.feature_scatter_plot.map {meta, feat_scatter_plot -> [feat_scatter_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.nfeature_plot.map {meta, nfeature_plot -> [nfeature_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.ncount_plot.map {meta, ncount_plot -> [ncount_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.cell_shape_plot.map{meta, cell_shape_plot -> [cell_shape_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.cell_segm_plot.map {meta, cell_segm_plot -> [cell_segm_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.area_histogram_plot.map {meta, area_histogram_plot -> [area_histogram_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.area_box_plot.map {meta, area_box_plot -> [area_box_plot]}.ifEmpty([]) )
        .combine ( SPATIAL_QC_POSTFILTER.out.area_overlapping_histogram_plot.map {meta, area_histogram_plot -> [area_histogram_plot]}.ifEmpty([]) )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'umap' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'split_cluster' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'dot_plot' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'vln_plot' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([file('NO_VLH')])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'umap' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'split_cluster' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'dot_plot' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'vln_plot' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'Harmony' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([file('NO_VAH')])
        )        
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'umap' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'split_cluster' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'dot_plot' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'vln_plot' &&
                    meta.normalization == 'log_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([file('NO_VLB')])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'umap' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'split_cluster' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'dot_plot' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([])
        )
        .combine (
            COMPILE_IMAGES_TO_VIDEO.out.image_video
                .filter { meta, video -> 
                    meta.plot_type == 'vln_plot' &&
                    meta.normalization == 'area_norm' &&
                    meta.clustering_method == 'BANKSY' 
                }
                .map {meta, video -> [video]}
                .ifEmpty([file('NO_VAB')])
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
