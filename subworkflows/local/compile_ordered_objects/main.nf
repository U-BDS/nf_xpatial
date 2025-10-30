#!/usr/bin/env nextflow
include { COMPILE_OBJECTS                          } from '../../../modules/local/compile_objects'

workflow COMPILE_ORDERED_OBJECTS {
    take:
        ch_xenium_obj // channel: xenium_objs

    main:
        ch_versions = Channel.empty()

        //
        // MODULE: Compile objects into a list
        //

        ch_xenium_obj
                .map { meta, xenium_obj -> [xenium_obj] }
                .collect()
                .map { xenium_obj_list -> 
                    [xenium_obj_list.sort { a, b ->
                        def filenameA = file(a).name
                        def filenameB = file(b).name
                        return filenameA.compareTo(filenameB)}]
                }

        COMPILE_OBJECTS (
            ch_xenium_obj
                .map{
                    meta, xenium_obj -> [xenium_obj]
                }
                .collect()
                .map{
                    [ [ 'id': 'compiled' ], it ]
                }
        )

    emit:
        compiled_obj = COMPILE_OBJECTS.out.compiled_obj
        versions     = ch_versions

}
