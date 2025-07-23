process MERGE_CLUSTER_TSV {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.1' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.1' 
        }"

    input:
    tuple val(meta), path(cluster_tsv)

    output:
    tuple val(meta), path("*.csv"), emit: merged_cluster_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ""
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)

    tsv_file_list <- str_split_1("${cluster_tsv}", "\\\\s")

    dfs <- list()

    for (tsv_file in tsv_file_list) {
        # Read the file
        df <- read.table(file = tsv_file, sep = '\\t', header = TRUE)

        # Use a unique identifier from the file path for the column name
        # path_parts <- strsplit(tsv_file, "/")[[1]]
        # unique_id <- path_parts[length(path_parts) - 1]  # Adjust based on the path structure

        # Rename the second column with the unique identifier
        # colnames(df)[2] <- unique_id

        # Append to the list of DataFrames
        dfs <- append(dfs, list(df))
    }

    # Print the number of DataFrames
    print(length(dfs))

    # Merge all DataFrames on the 'Index' column
    compiled_df <- Reduce(function(x, y) full_join(x, y, by = "Index"), dfs)

    #compiled_df <- compiled_df %>%
    #rename_with(
    #    ~ paste0("clust_", .),
    #    matches("^(AreaNorm|Xenium|spe)_BSKY_")
    #)

    # Save the result to a CSV file if needed
    write.csv(compiled_df, "${prefix}_compiled_clusts_results.csv")
    """
}
