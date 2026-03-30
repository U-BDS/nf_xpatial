process COMPILE_OBJECTS {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
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
    file_list <- sort(file_list)
    
    obj_list = list()

    for (i in 1:length(file_list)) {
        xenium.obj <- readRDS(file = file_list[[i]])

        obj_list[[i]] <- xenium.obj
    }

    print(obj_list)

    saveRDS(
        object = obj_list,
        file = "${prefix}_compiled.rds"
    )
    """
}
