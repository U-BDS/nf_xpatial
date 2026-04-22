# U-BDS/nf_xpatial: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### Output Files

-   `<sample_name>/`

    -   `filtered/`

        -   `<sample_name>.filtered.csv` : This contains the amount of cells filtered based on the different criteria implemented in the pipeline

        -   `<sample_name>.filtered.rds` : This is the post-filtered Seurat object for this sample

    -   `manual_annotations/`

        -   `<sample_name>_annotated.rds` : This contains the Seurat object that contains any manual annotations from the Xenium explorer. This file appears if you enter region annotations or regions to remove from the tissue.

    -   `normalized/`

        -   `area/`

            -   `<sample_name>_norm_area.rds` : This is the post-area normalized Seurat object.

        -   `log/`

            -   `<sample_name>_norm_log.rds` : This is the post-log normalized Seurat object.

    -   `raw_data/`

        -   `<sample_name>.rds` : This is the raw Seurat object, created directly from the Xenium outputs

-   `compiled/`

    -   `BANKSY/` or `Harmony/` : The structure is the same for these two subdirectories. The images act as a way to see the effects of the different clustering algorithms and parameters to determine the best tool and parameter choice to use for downstream analysis.

        -   `area_norm/` or `log_norm/` : The structure is the same for these two subdirectories. The code is run separately for each normalization method.

            -   `cluster_csvs/`

                -   `*.csv` : Each csv file contains the clustering information for a single clustering parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) that can be added to a Seurat object as metadata.

            -   `qc/`

                -   `dot_plot/`

                    -   `*.pdf` : This contains the dot plots for every marker provided in the marker list. Each group is present on a separate page. Each parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) is a separate file.

                -   `split_cluster_plots/`

                    -   `*.png` : This is a faceted plot with clusters split out by sample, projected onto the Spatial Dim plot, and showing the proportions for each sample. Each parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) is a separate file.

                -   `umap_dim_plot/`

                    -   `*.png` : This is the umap plot. Each parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) is a separate file.

                -   `vln_plot/`

                    -   `*.pdf` : This contains the vln plots for every marker provided in the full marker list. Each group is present on a separate page. Each parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) is a separate file.

            -   `reduced_dim_csvs/`

                -   `*.csv` : Each csv file contains the dimension reduction information for a single parameter combination (dimension and resolution for `Harmony/` or lambda, k_geom, dimension, and resolution for `BANKSY/`) that can be added to a Seurat object as a new Dimension Reduction. To add a new Dimension Reduction, you'll need to add the loadings, embeddings, and stdev files.

            -   `*.rds` : Each object contains the clustering information and dimension reductions for a single parameter combination. The parameter combination are present as part of the file name.

    -   `filtered/` or `raw_data/` : These directories contain the same structure. We generate all the plots to show the differences (if there any) between the same metrics before (`raw_data/`) and after filtering (`filtered/`).

        -   `qc/`

            -   `cell_area_qc/`

                -   `compiled_filtered_box_plot.png` : This shows cell ares for each cell into a box plot split out by sample.

                -   `compiled_filtered_overlapping_histogram_plot.png` : This compiles the cell areas for cell into a single histogram plot, with each sample overlapping one another.

                -   `compiled_histogram_plot.png` : This is a tiled histogram showing cell area for each cell, each facet of the plot is a single sample.

            -   `cell_shape_qc/`

                -   `compiled.cell_segmentation_proportion_plot.png` : This plot shows the proportion of the cell segmentation method for each cell that comprise each sample.

                -   `compiled.cell_shape_proportion_plot.png` : This plot shows the proportion of cell shape that comprise each sample.

            -   `compiled_dim_plot.png` : Dim plots for the samples, coloring based on Biological Group.

            -   `compiled_feature_scatter_plot.png` : These show nCount vs nFeature for each sample

            -   `compiled_vln_plot.png` : This is a tiled plot, showing nCount and nFeature for each sample to determine reasonable cutoffs.

            -   `compiled.qc_ncount_image_feature_plot.png` : Shows nCount distribution across the tissue for all samples.

            -   `compiled.qc_nfeature_image_feature_plot.png` : Shows nFeature distribution across the tissue for all samples.

    -   `manual_annotations/`

        -   `qc/`

            -   `compiled_annotated_dim_plot.png` : This is a dim plot for samples that had manual annotations. In this plot, we would show the anatomical/pathological annotations mapped onto the Seurat objects. This is also generated if sections of a tissue are removed via labelling them to be removed.

    -   `compiled_area_norm_all_clusters.rds` : This is one of the final objects output by the pipeline. This contains all the Seurat and BANKSY clustering parameter combinations as metadata and all dimensions for the area normalized + integrated seurat object

    -   `compiled_log_norm_all_clusters.rds` : This is one of the final objects output by the pipeline. This contains all the Seurat and BANKSY clustering parameter combinations as metadata and all dimensions for the log normalized + integrated seurat object

