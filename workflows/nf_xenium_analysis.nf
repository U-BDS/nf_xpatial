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
include { CREATE_XENIUM_OBJ                         } from '../modules/local/create_xenium_object'
include { MAKE_IMAGE_DIM_PLOT as RAW_IMAGE_DIM_PLOT } from '../modules/local/make_image_dim_plot'
include { COMPILE_OBJECTS                           } from '../modules/local/compile_objects'
include { ADD_MANUAL_ANNOTATIONS                    } from '../modules/local/add_manual_annotations'

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
    CREATE_XENIUM_OBJ(
        ch_samplesheet
            .map{
                meta, xenium_input, metadata, manual_annotation ->
                [meta, xenium_input, metadata]
            }
    )
    ch_versions = ch_versions.mix(CREATE_XENIUM_OBJ.out.versions)

    // START ADD_MANUAL_ANNOTATION SUBWORKFLOW
    
    // Separate the samples that have manual annotations
    ch_samplesheet
        .join(CREATE_XENIUM_OBJ.out.xenium_obj)
        .map{
            meta, xenium_input, metadata, manual_annotation, xenium_rds ->
                [meta, xenium_rds, manual_annotation]
        }
        .view()
        .branch{
            meta, xenium_rds, manual_annotation ->
                with_annotation: manual_annotation
                no_annotation: true
        }
        .set{ ch_sep_objects }
    
    //
    // MODULE: Add manual annotations where possible
    //

    ADD_MANUAL_ANNOTATIONS(
        ch_sep_objects.with_annotation
    )
    ch_annotated_xenium_obj = ADD_MANUAL_ANNOTATIONS.out.annotated_xenium_obj

    // END ADD_MANUAL_ANNOTATION SUBWORFLOW

    //
    // MODULE: Compile objects into a list
    //
    COMPILE_OBJECTS(
        CREATE_XENIUM_OBJ.out.xenium_obj
            .map{
                meta, xenium_obj -> [xenium_obj]
            }
            .collect()
            .map{
                [ [ 'id': 'compiled_RAW' ], it ]
            }
    )

    //
    // MODULE: Create an initial Image Dim Plot
    //
    RAW_IMAGE_DIM_PLOT(
        COMPILE_OBJECTS.out.compiled_obj
    )
    ch_versions = ch_versions.mix(RAW_IMAGE_DIM_PLOT.out.versions)


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
