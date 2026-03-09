#!/usr/bin/env nextflow

include { FIND_VARIABLE_FEATURES } from '../../../modules/local/find_variable_features'
include { COMPILE_LISTS          } from '../../../modules/local/compile_lists'

workflow GET_VARIABLE_FEATURES {
    take:
        ch_xenium_obj           // channel: xenium objects

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Find variable features
        //
        FIND_VARIABLE_FEATURES ( 
            ch_xenium_obj 
        )

        //
        // MODULE: Compile variable feature lists
        //
        COMPILE_LISTS (
            FIND_VARIABLE_FEATURES.out.variable_feature_list
                .map { meta, gene_list -> [gene_list] }
                .collect()
                .map { 
                    gene_list -> 
                        def new_meta = ['id': 'genes']
                    [new_meta, gene_list]
                }
        )

    emit:
        versions      = ch_versions
        vf_xenium_obj = FIND_VARIABLE_FEATURES.out.variable_features_xenium_obj
        gene_list     = COMPILE_LISTS.out.compiled_list.map { meta, gene_list -> [gene_list] }


}
