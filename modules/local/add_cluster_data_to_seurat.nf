process ADD_CLUSTER_DATA_TO_SEURAT {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(xenium_object), val(cluster_metadata_csv), val(embeddings_csv), val(loadings_csv), val(stdev_csv)

    output:
    tuple val(meta), path("*.rds"), emit: all_cluster_xenium_obj
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = "--assay ${meta.assay}"

    """
    add_cluster_data_to_seurat.R \\
        $args \\
        $assay_flag \\
        --input "$xenium_object" \\
        --clusters "${cluster_metadata_csv}" \\
        --embeddings "${embeddings_csv.join(',')}" \\
        --loadings "${loadings_csv.join(',')}" \\
        --stdev "${stdev_csv.join(',')}" \\
        --outfile ${prefix}_all_clusters.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
