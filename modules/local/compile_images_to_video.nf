process COMPILE_IMAGES_TO_VIDEO {
    tag "$meta.id"
    label 'process_high'

    container "${ 
        (workflow.containerEngine == 'singularity') &&
            (!task.ext.singularity_pull_docker_container) ?
            'docker://uabbds/nf_xenium_analysis:0.0.4' :
            'docker.io/uabbds/nf_xenium_analysis:0.0.4' 
        }"

    input:
    tuple val(meta), path(image_list)

    output:
    tuple val(meta), path("*.mp4"), emit: image_video

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(av)
    library(pdftools)
    library(stringr)

    file_list <- str_split_1("${image_list}", "\\\\s")

    plots_png <- file_list

    # Convert to png
    if (any(endsWith(file_list, ".pdf"))) {
        plots_png <- lapply(
            file_list,
            function(x) {
                pdf_convert(x, format = "png", pages = NULL, dpi = 300, filenames = gsub("pdf", "png", x))
            }
        )
    }

    av_encode_video(
        unlist(plots_png),
        framerate = 0.5,
        output = "${prefix}.mp4"
    )

    """
}
