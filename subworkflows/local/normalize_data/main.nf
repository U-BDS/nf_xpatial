#!/usr/bin/env nextflow

include { NORMALIZE_LOG                                 } from '../../../modules/local/normalize_log'
include { COMPILE_OBJECTS as COMPILE_LOG_OBJECTS        } from '../../../modules/local/compile_objects'
//include { NORMALIZE_SCTRANSFORM                         } from '../../../modules/local/normalize_sctransform'
//include { COMPILE_OBJECTS as COMPILE_SCT_OBJECTS        } from '../../../modules/local/compile_objects'
include { NORMALIZE_AREA                                } from '../../../modules/local/normalize_area'
include { COMPILE_OBJECTS as COMPILE_AREA_OBJECTS       } from '../../../modules/local/compile_objects'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NORM_NCOUNT   } from '../../../modules/local/qc_image_feature_plot'
include { QC_IMAGE_FEATURE_PLOT as QC_IFP_NORM_NFEATURE } from '../../../modules/local/qc_image_feature_plot'

workflow NORMALIZE_DATA {
    take:
        ch_xenium_obj         // channel: xenium objects
        normalization_methods // list: list of user-selected normalization methods

    main:
        ch_versions = Channel.empty()
        ch_normalized_objects = Channel.empty()
        ch_compiled_norm_objects = Channel.empty()

        //
        // MODULE: Log Normalization
        //
        if ('log' in normalization_methods){
            NORMALIZE_LOG (
                ch_xenium_obj
            )

            ch_versions = ch_versions.mix(NORMALIZE_LOG.out.versions)
            ch_log_norm_objs = NORMALIZE_LOG.out.normalized_xenium_obj
                .map {
                    meta, xenium_obj ->
                        def new_meta = meta + [normalization: 'log_norm']
                    [new_meta, xenium_obj]
                }

            ch_normalized_objects = ch_normalized_objects.mix (ch_log_norm_objs)

            COMPILE_LOG_OBJECTS (
                ch_log_norm_objs
                    .map{
                        meta, xenium_obj -> [xenium_obj]
                    }
                    .collect()
                    .map{
                        [ [ 'id': 'compiled_log_norm', 'normalization': 'log_norm' ], it ]
                    }
            )
            ch_compiled_norm_objects = ch_compiled_norm_objects.mix(COMPILE_LOG_OBJECTS.out.compiled_obj)
        }

        //
        // MODULE: SCTransform Normalization
        //
        //if ('sctransform' in normalization_methods) {
        //    NORMALIZE_SCTRANSFORM (
        //        ch_xenium_obj
        //    )
        //    
        //    ch_versions = ch_versions.mix(NORMALIZE_SCTRANSFORM.out.versions)
        //    ch_sct_norm_objs = NORMALIZE_SCTRANSFORM.out.normalized_xenium_obj
        //        .map {
        //            meta, xenium_obj ->
        //                def new_meta = meta + [normalization: 'sctransform']
        //            [new_meta, xenium_obj]
        //        }

        //    ch_normalized_objects = ch_normalized_objects.mix (ch_sct_norm_objs)

        //    COMPILE_SCT_OBJECTS (
        //        ch_sct_norm_objs
        //            .map{
        //                meta, xenium_obj -> [xenium_obj]
        //            }
        //            .collect()
        //            .map{
        //                [ [ 'id': 'compiled_sct_norm' ], it ]
        //            }
        //    )
        //    ch_compiled_norm_objects = ch_compiled_norm_objects.mix(COMPILE_SCT_OBJECTS.out.compiled_obj)
        //}

        //
        // MODULE: Area Normalization
        //
        if ('area' in normalization_methods) {
            NORMALIZE_AREA (
                ch_xenium_obj
            )

            ch_versions = ch_versions.mix(NORMALIZE_AREA.out.versions)
            ch_area_norm_objs = NORMALIZE_AREA.out.normalized_xenium_obj
                .map {
                    meta, xenium_obj ->
                        def new_meta = meta + [normalization: 'area_norm']
                    [new_meta, xenium_obj]
                }

            ch_normalized_objects = ch_normalized_objects.mix (ch_area_norm_objs)

            COMPILE_AREA_OBJECTS (
                ch_area_norm_objs
                    .map{
                        meta, xenium_obj -> [xenium_obj]
                    }
                    .collect()
                    .map{
                        [ [ 'id': 'compiled_area_norm', 'normalization': 'area_norm'], it ]
                    }
            )
            ch_compiled_norm_objects = ch_compiled_norm_objects.mix(COMPILE_AREA_OBJECTS.out.compiled_obj)
        }

        ch_ifp_nfeature_in = ch_compiled_norm_objects
        ch_ifp_ncount_in = ch_compiled_norm_objects

        //
        // MODULE: nCount Feature Plot
        //
        QC_IFP_NORM_NFEATURE(
            ch_ifp_nfeature_in
        )
        ch_versions = ch_versions.mix(QC_IFP_NORM_NFEATURE.out.versions)


        //
        // MODULE: nFeature Feature Plot
        //
        QC_IFP_NORM_NCOUNT(
            ch_ifp_ncount_in
        )
        ch_versions = ch_versions.mix(QC_IFP_NORM_NCOUNT.out.versions)


    emit:
        normalized_objects          = ch_normalized_objects
        compiled_norm_objects       = ch_compiled_norm_objects
        image_feature_plot_ncount   = QC_IFP_NORM_NCOUNT.out.image_feature_plot
        image_feature_plot_nfeature = QC_IFP_NORM_NFEATURE.out.image_feature_plot
        versions                    = ch_versions


}
