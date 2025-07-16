process EXTRACT_REDUCED_DIMS {
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
    tuple val(meta), path(spe_obj)

    output:
    tuple val(meta), path("*.tsv"), emit: reduced_dims_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assay_flag = meta.normalization == 'area_norm' ? '--assay AreaNorm' : '--assay Xenium'

    """
    extract_reduced_dims.R \\
        $args \\
        $assay_flag \\
        --input "$spe_obj" \\
        --outfile "${prefix}_reduced_dims.tsv"
    """
}
