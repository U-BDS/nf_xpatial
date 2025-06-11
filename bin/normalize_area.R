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
        c("-a", "--cell_area_norm_assay"),
        type="character",
        default="AreaNorm",
        help="The name of the assay in the metadata that will contain the cell area norm information."),
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
expr_mtx <- GetAssayData(object = xenium_obj, assay = opt$assay, layer = opt$layer)
counts <- t(expr_mtx)

# Confrim cell order matches between metadata and counts mtx
if (!all(rownames(counts) == rownames(xenium_obj@meta.data))) {
    stop("Cell Ids in metadata do not match those in counts")
}

# Normalize counts by scaling factor
norm_counts <- t(counts / scaling_factor)

# Apply log transformation
log_norm_counts <- log1p(norm_counts)

# Store the data on the xenium_obj
xenium_obj[[opt$cell_area_norm_assay]] <- CreateAssayObject(counts = norm_counts)
LayerData(xenium_obj, assay = opt$cell_area_norm_assay, layer = "data") <- log_norm_counts

# Set default assay
DefaultAssay(xenium_obj) <- opt$cell_area_norm_assay

xenium_obj <- FindVariableFeatures(xenium_obj)

xenium_obj <- ScaleData(xenium_obj)

DefaultAssay(xenium_obj) <- "Xenium"

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
