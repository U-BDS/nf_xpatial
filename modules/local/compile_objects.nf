process COMPILE_OBJECTS {
    tag "$meta.id"
    label 'process_low'

    //container "nf_xenium_analysis_0.0.1.sif"
    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.1' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.1' 
        }"

    input:
    tuple val(meta), path(xenium_list)

    output:
    tuple val(meta), path("*.rds"), emit: compiled_obj

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(Seurat)
    library(stringr)

    file_list <- str_split_1("${xenium_list}", "\\\\s")
    
    obj_list = list()

    for (i in 1:length(file_list)) {
        xenium.obj <- readRDS(file = file_list[[i]])

        obj_list[[i]] <- xenium.obj
    }

    print(obj_list)

    saveRDS(
        object = obj_list,
        file = "${prefix}_merged.rds"
    )

    """
}
