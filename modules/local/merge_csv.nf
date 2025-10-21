process MERGE_CSV {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(cluster_csv)

    output:
    tuple val(meta), path("*.csv"), emit: merged_cluster_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ""
    def prefix     = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)

    csv_file_list <- str_split_1("${cluster_csv}", "\\\\s")

    dfs <- list()

    for (csv_file in csv_file_list) {
        # Read the file
        df <- read.table(file = csv_file, sep = ',', header = TRUE)

        # Append to the list of DataFrames
        dfs <- append(dfs, list(df))
    }

    # Print the number of DataFrames
    print(length(dfs))

    # Merge all DataFrames on the 'Index' column
    compiled_df <- Reduce(function(x, y) full_join(x, y, by = "Index"), dfs)

    # Save the result to a CSV file if needed
    write.csv(compiled_df, "${prefix}_merged_results.csv", quote = FALSE, row.names = FALSE)
    """
}
