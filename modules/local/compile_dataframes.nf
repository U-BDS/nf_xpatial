process COMPILE_DATAFRAMES {
    tag "$meta.id"
    label 'process_low'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(df_list)

    output:
    tuple val(meta), path("*.csv"), emit: compiled_df

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(stringr)

    file_list <- str_split_1("${df_list}", "\\\\s")
    
    df_list = list()

    for (i in 1:length(file_list)) {
        df_list[[i]] <- read.csv(
            file = file_list[[i]],
            header = TRUE,
            row.names = 1
        )
    }

    combined_df <- do.call(rbind, df_list)

    write.csv(
        combined_df,
        file = "${prefix}_compiled.csv",
        row.names = TRUE
    )

    """
}
