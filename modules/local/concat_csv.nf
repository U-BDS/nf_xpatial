process CONCAT_CSV {
    tag "$meta.id"
    label 'process_medium'

    //container "nf_xenium_analysis_0.0.1.sif"
    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.1' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.1' 
        }"

    input:
    tuple val(meta), path(csv_list)

    output:
    tuple val(meta), path("*.csv"), emit: concat_csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env bash

    out_file="${prefix}_merged.csv"

    for csv_file in ${csv_list}
    do
        # print header to file
        if [[ ! -f "\${out_file}" ]]; then
            head -n1 "\${csv_file}" > "\${out_file}"
        fi

        tail -n +2 "\${csv_file}" >> "\${out_file}"
        
    done

    """
}
