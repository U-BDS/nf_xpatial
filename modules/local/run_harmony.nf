process RUN_HARMONY {
    tag "$meta.id"
    label 'process_medium'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'library://atrull314/uabbds/nf_xpatial:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(xenium_object)

    output:
    tuple val(meta), path("*.rds"), emit: integrated_xenium_obj
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ""
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def assay_flag = "--assay ${meta.assay}"

    """
    run_harmony.R \\
        $args \\
        $assay_flag \\
        --reduction_name "pca_${meta.normalization}" \\
        --input "$xenium_object" \\
        --outfile ${prefix}_harmony.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
