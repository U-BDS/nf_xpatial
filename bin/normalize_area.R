#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(arrow)      # For reading/writing Parquet files
library(dplyr)      # For data manipulation
library(jsonlite)   # For working with JSON data
library(knitr)      # For generating reports
library(optparse)   # For parsing commandline arguments
library(progressr)  # For progress bars
library(purrr)      # Functional programming tools
library(Seurat)     # Main analysis package

# Plotting
library(patchwork)  # For combining plots

# Set options
options(future.globals.maxSize = 8192 * 1024^2)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
  make_option(
    c("-i", "--input"),
    type="character",
    default=NULL,
    metavar="path",
    help="The xenium results to be analyzed"),
  make_option(
    c("-a", "--assay"),
    type="character",
    default="Xenium",
    help="The assay to operate on"),
  make_option(
    c("-l", "--counts_layer"),
    type="character",
    default="counts",
    help="The layer in the assay with raw counts"),
  make_option(
    c("-c", "--cell_area_col"),
    type="character",
    default="Cell_Area",
    help="The name of the column in the metadata containing the cell areas"),
  make_option(
    c("-n", "--cell_area_norm_assay"),
    type="character",
    default="AreaNorm",
    help="The name of the assay in the metadata that will contain the cell area norm information."),
  make_option(
    c("-r", "--vars_to_regress"),
    type="character",
    default=NULL,
    help="The variables to regress out during scaling"),
  make_option(
    c("-f", "--nfeatures"),
    type="integer",
    default=2000,
    help="Number of features to select as top variable features"),    
  make_option(
    c("-o", "--outfile"),
    type="character",
    default="normalized_xenium_obj.rds",
    metavar="path",
    help="The output name for the seurat object"))

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Check that the cell area column exists
if ( !opt$cell_area_col %in% colnames(xenium_obj@meta.data) ) {
    stop("Unable to find the column containing cell area", call. = FALSE)
}

#################
### NORMALIZE ###
#################

# Compute the scaling factor for area normalization
scaling_factor <- xenium_obj@meta.data[[opt$cell_area_col]] / median(xenium_obj@meta.data[[opt$cell_area_col]])

# Extract the expression matrix
expr_mtx <- GetAssayData(object = xenium_obj, assay = opt$assay, layer = opt$counts_layer)
counts <- t(expr_mtx)

# Confirm cell order matches between metadata and counts mtx
if (!all(rownames(counts) == rownames(xenium_obj@meta.data))) {
    stop("Cell Ids in metadata do not match those in counts")
}

# Normalize counts by scaling factor
norm_counts <- t(counts / scaling_factor)

# Apply log transformation
log_norm_counts <- log1p(norm_counts)

# Store the data on the xenium_obj 
# `count` == area normalized counts, `data` == log transformed area normalized counts
xenium_obj[[opt$cell_area_norm_assay]] <- CreateAssayObject(counts = norm_counts)
LayerData(xenium_obj, assay = opt$cell_area_norm_assay, layer = "data") <- log_norm_counts

# Set default assay
DefaultAssay(xenium_obj) <- opt$cell_area_norm_assay

# warning when nfeatures >= than total number of features available
num_features <- opt$nfeatures

if (opt$nfeatures > nrow(xenium_obj[[opt$assay]]$counts)) {
    warning(paste0(
        "The number of total features available in the current assay (",
        nrow(xenium_obj[[opt$assay]]$counts),
        ") is less than ", 
        opt$nfeatures,
        "\nFindVariableFeatures will results in using all available feaures."
        )
    )

    num_features <- nrow(xenium_obj[[opt$assay]]$counts)
}

xenium_obj <- FindVariableFeatures(xenium_obj,
                                   nfeatures = num_features)

xenium_obj <- ScaleData(
  xenium_obj,
  vars.to.regress = opt$vars_to_regress
)

#################
### SAVE DATA ###
#################

saveRDS(
    object = xenium_obj,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
