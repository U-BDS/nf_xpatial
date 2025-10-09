process CLASSIFY_CELL_SHAPE{
    tag "$meta.id"
    label 'process_low'

    //container "nf_xenium_analysis_0.0.2.sif"
    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.2' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.2' 
        }"

    input:
    tuple val(meta), path(xenium_obj)

    output:
    tuple val(meta), path('*.rds'), emit: cell_shape_xenium_obj
    tuple val(meta), path('*.csv'), emit: cell_shape_csv
    path 'versions.yml'           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    classify_cell_shape.R \\
        $args \\
        --input $xenium_obj \\
        --outfile ${prefix}_cell_shape.rds \\
        --shape_file ${prefix}_cell_shape.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        r-seurat: \$(Rscript -e "library(Seurat); cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
