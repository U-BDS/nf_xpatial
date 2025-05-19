process ADD_MANUAL_ANNOTATIONS{
    tag "$meta.id"
    label 'process_low'

    container "nf_xenium_analysis_0.0.1.sif"

    input:
    tuple val(meta), path(xenium_object), path(manual_annotations)

    output:
    tuple val(meta), path("*.rds"), emit: annotated_xenium_obj
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    add_manual_annotation.R \\
        $args \\
        --input "$xenium_object" \\
        --manual_annotation $manual_annotations \\
        --outfile ${prefix}_annotated.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
