![nf_xpatial](./assets/nf_xpatial_logo_transparent.png)

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

- [nf\_xpatial](#u-bdsnf_xpatial)
  - [Introduction](#introduction)
  - [Pipeline Summary](#pipeline-summary)
  - [Usage](#usage)
  - [Outputs](#outputs)
  - [Detailed outputs](#detailed-outputs)
  - [Notes](#notes)
  - [Credits](#credits)
  - [Contributions and Support](#contributions-and-support)
  - [Citations](#citations)


## Introduction

**nf_xpatial** is a best-practices bioinformatics pipeline written in Nextflow that can be used to perform tertiary analysis on 10x Xenium data. It uses the output directories produced by the Xenium Onboard Analysis (XOA) instrument as input and performs quality control, filtering, normalization, and clustering, and generates configurable figures that can be reviewed individually or in the final summary report.

Notably, for cases where the raw data no longer follows the Xenium Onboard Analysis outputs (e.g.: sample was re-segmented with a third parity tool) nf_xpatial also accepts Seurat objects as input (per sample Seurat object without any data processed).

## Pipeline Summary

![pipeline diagram](./assets/nf_xpatial_metro.png)

1. Create [Seurat](https://github.com/satijalab/seurat) object(s) from Xenium output
2. Generate QC images for raw data
   1. Cell Area QC (`Area Box Plot`, `Area Histogram Plot`, `Overlapping Histogram Plot`)
   2. Cell Shape QC (`Cell Segmentation Proportion Plot`, `Cell Shape Proportion Plot`)
   3. General QC (`Image Dim Plot`, `nFeature/nCount Violin Plot`, `nFeature/nCount Feature Scatter Plot`, `nFeature Dim Plot`, `nCount Dim Plot`)
3. Filter the Seurat object
4. Generate QC images for post-filetered data
   1. Cell Area QC (`Area Box Plot`, `Area Histogram Plot`, `Overlapping Histogram Plot`)
   2. Cell Shape QC (`Cell Segmentation Proportion Plot`, `Cell Shape Proportion Plot`)
   3. General QC (`Image Dim Plot`, `nFeature/nCount Violin Plot`, `nFeature/nCount Feature Scatter Plot`, `nFeature Dim Plot`, `nCount Dim Plot`)
5. Normalize the data: choose between [`area normalization`](./docs/usage.md#area-normalize), [`log normalization`](./docs/usage.md#log-normalize), or execute both
6. Gene Pair QC (`Barnyard Plot`, `Heatmap Plot`)
7. Merge normalized Seurat objects
8. Integrate the data (with [Harmony](https://github.com/immunogenomics/harmony))
8. Perform Seurat clustering for single-cell clustering
   1. Scale data
   2. Run PCA
   3. Run Harmony
   4. Run UMAP
   5. Find Neighbors
   6. Find Clusters
9. Perform [BANKSY](https://github.com/prabhakarlab/Banksy) clustering for single-cell and/or spatial-domain clustering (**Note**: this can be executed with BANKSY or with the BANKSY Seurat Wrapper)
   1. Convert to Spatial Experiment
   2. Stagger Spatial Coordinates
   3. Compute BANKSY Matrix
   4. Compute BANKSY PCA 
   5. Run Harmony BANKSY
   6. Run BANKSY UMAP
10. Merge BANKSY and Seurat clustered objects into a single object
11. Generate Cluster QC images (**This is done for all parameter combinations**) (`UMAP Dim Plot`, `Split Cluster Plot`, `Marker Violin Plot`, `Marker Dot Plot`)
12. Generate summary report

## Usage

### Create a metadata.csv (required)

First, prepare a csv file containing metadata for the samples to be analyzed. A user can choose to create separate metadata csv's for each sample or create a single metadata csv that contains information for all samples. The only required columns in this file are `SampleID` and `BiologicalGroup`, however additional columns can be added which will be stored in the Seurat object.

`metadata.csv`

```csv
SampleID,BiologicalGroup
XNM001,Control
XNM002,Treatment
XNM003,Control
XNM004,Treatment
```

### Create Tissue annotation files (optional)

If you have any tissue annotations, i.e. regions that you have drawn and labelled using the Xenium Explorer, you are able to add these onto the seurat object. Once the annotations are exported outside of Xenium Explorer, these will need to be reformatted so that all the annotations are in a tab-delimited file with the columns `Cell_ID` and `Tissue_annotation`. To assist with this step, we provide a script in this repository (`bin/gather_xenium_explorer_annotations.sh`) that can be used to process the exports from Xenium Explorer into the format needed by this pipeline.

Additionally, this step can also be used to remove parts of a slide by labelling the region you wish to remove as `REMOVE` in Xenium Explorer. The most common use cases for this are to remove parts of sample that has folded over on itself or to remove regions that are from a different sample (NOTE: The pipeline does not currently have a way to add these regions back to the sample it belongs to).

An example file format for the cases described above is presented below:

Case where specific cells map to regions in slide called `a`, `b`, `c`

```csv
Cell_ID	Tissue_annotation
efhphlac-1	a
fpldnmpm-1	c
bmhpfjfb-1	c
eonmgbhj-1	b
gmjbldbh-1	c
cjfbfjmn-1	b
```

Case where cells below will be removed from any downstream processing

```csv
Cell_ID	Tissue_annotation
mbohedjb-1	REMOVE
mbohdkcj-1	REMOVE
bieojgni-1	REMOVE
mbpkcimm-1	REMOVE
mbpkhbip-1	REMOVE
mbohjhae-1	REMOVE
```

### Create a samplesheet.csv (required)

Finally, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,xenium,metadata,manual_annotation
XNM001,/path/to/XNM001_xenium_output,/path/to/xenium_metadata.csv,/path/to/XNM001_manual_annotation.csv
XNM002,/path/to/XNM002_xenium_output,
XNM003,/path/to/XNM003_xenium_output,
XNM004,/path/to/XNM004_xenium_output,/path/to/xenium_metadata.csv,/path/to/XNM004_manual_annotation.csv
```

Each row represents a directory produced by the Xenium Onboard Analysis.

### Create a gene marker list (recommended)

A marker list is recommended to support the evaluation of the single-cell and spatial clusters. This list should be a 2 column csv file such as the example below and provided to the `--marker_gene_list` option of nf_xpatial:

```csv
group,gene
Celltype Marker,Gad1
Celltype Marker,Drd1
Celltype Marker,Drd2
Celltype Marker,Slc17a6
Celltype Marker,Slc17a7
Celltype Marker,Aqp4
Celltype Marker,Pdgfra
```

Although the workflow supports multiple categories in the `group` column, only the first category in the csv is displayed in the summary report for simplicity (the example above uses a single category for all cell type markers). While only the first group is present in the report, in the case multiple groups are given, all groups are shown in the figures generated (dot plots and violin plots will group genes based on the 'group' column and will print each 'group' on a separate page in their pdf output). By default they will separate groups if they are over 50 genes long, and will spread the group out over multiple pages if it is. This can be configured via custom config by modifying the `--max_genes_per_group <INT>` parameter for the `MARKER_DOT_PLOT` and `MARKER_VLN_PLOT` processes.

### Execute the workflow

Now, you can run the pipeline using (in the example below both `area` and `log` normalization methods are enabled, but a user can choose only one of them):

```bash
nextflow run U-BDS/nf_xpatial \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --normalization_method "area,log" \
   --dim_Seurat "25,30" \
   --res_Seurat "0.4,0.5,0.6,0.7" \
   --lambda_BANKSY "0.2,0.8" \
   --k_geom_BANKSY "15,30" \
   --nPCs_BANKSY "20,30" \
   --res_BANKSY "0.4,0.6,0.8,1.0" \
   --outdir <OUTDIR>
```

For more details on enabling additional parameters, or usage please refer the advanced [usage documentation](/docs/usage.md).

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;

## Outputs

![Sample_param_compilation](/docs/images/nf_xpatial_param_comp.gif)
(Sample data from Jeremy Day, Jamie Peters and Jasper Heinsbroek)

`nf_xpatial` produces a number of files and figures that can be used to review the quality of the data and refine clustering. However, the main output of this pipeline are `.rds` objects that contain all clustering results into a central object. An `.rds` object is created for each normalization method specified by the `--normalization_method` pipeline parameter.

Because these objects contain all clustering results and all dimension reductions they can be quite large, making it prudent to filter these objects to a single (or selection) of parameter combinations. In order to do that, it's important to note how the data is stored on each object:

1. Each normalization stores its result in a specific assay, `log_norm` stores its data in the `Xenium` assay, while `area_norm` stores its data in the `AreaNorm` assay. 
2. The clustering calls are stored in the objects metadata with the following format: `clust_[method]_[clustering parameters]`. `method` may be one of 1. `SEU` (Seurat clustering) 2. `BSKY` (BANKSY clustering) 3. `BSKYSEU`(BANSKY’s Seurat wrapper). `clustering parameters` may be `d` (PCA dimensions) , `r` (resolution), `l` (lambda), `k` (k_geom).
3. UMAP dimension reductions follow a similar format to clustering, specifically `[method]_[reduction]_[clustering parameters]`. `method` is the same as described above, `reduction` may be one of 1. `pca` 2. `harmony`, 3. `umap` and the `parameters` option match those described above, with the exception that reductions are calculated prior to clustering, so there are no resolution (`r`) parameters in their names.

### Detailed outputs

1. For in-depth descriptions and locations of the additional outputs within the results folder, refer to this [document](docs/output.md) (located at `docs/output.md`)

2. We provide a brief guide that details the naviation of the compiled objects which are produced by `nf_xpatial`. [This guide can be found here](/docs/nf_xpatial_navigation.md). To assist with filtering the object(s), we provided this [script](assets/filter_xenium_obj.R) (located at `assets/filter_xenium_obj.R`) to perform the filtering as well as listing some examples on how to use the provided script.

## Credits

U-BDS/nf_xpatial was originally written by [Luke Potter](https://github.com/LPotter21), [Nilesh Kumar](https://github.com/nilesh-iiita), [Austyn Trull](https://github.com/atrull314), [Lara Ianov](https://github.com/lianov).

We would also like to thank the following people and groups for their support, including financial support:

- Elizabeth Worthey
- Jeremy Day
- Jamie Peters
- Jasper Heinsbroek
- Frances Lund
- Funding: 
   - Health Services Foundation’s General Endowment Fund
   - University of Alabama at Birmingham Biological Data Science Core (U-BDS), RRID:SCR_021766, <https://github.com/U-BDS>
   - Civitan International Research Center
   - UAB Office of Research
   - 3P30CA013148-48S8
   - UM1TR004771
   - UAB MULTIPI8110 and Dr. Worthey's start-up funds

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use U-BDS/nf_xpatial for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

In addition, an extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
