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
        c("-a", "--manual_annotation"),
        type="character",
        default=NULL,
        metavar="path",
        help="The manual annotations file"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="xenium_obj.rds",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("--unknown"),
        type="character",
        default="UN",
        help="Annotation name to set unknown cells")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}

if (is.null(opt$manual_annotation)) {
    print_help(opt_parser)
    stop("Please provide the metadata file to add to the seurat object.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Read in manual_annotations
manual_annotation <- read.csv(opt$manual_annotation, sep = "\t", row.names = "Cell_ID")["Tissue_annotation"]

# Get list of cells from xenium
cell_ids <- Seurat::Cells(xenium_obj)

# Set missing cells to unknown value
missing_cells <- setdiff(cell_ids, rownames(manual_annotation))

if (length(missing_cells) > 0) {
    cat(
        length(missing_cells),
        "cells are missing in manual annotations. Setting their annotation to ", opt$unknown, ".\n",
        sep=""
    )

    missing_cell_df <- data.frame(
        Tissue_annotation = rep(opt$unknown, length(missing_cells)),
        row.names = missing_cells
    )
    manual_annotation <- rbind(manual_annotation, missing_cell_df)
} else {
    print("No cells are missing a manual annotation!\n")
}

# Reorder manual_annotations
manual_annotation <- manual_annotation[cell_ids, , drop = FALSE]

# Set Tissue annotations factors
manual_annotation$Tissue_annotation <- factor(manual_annotation$Tissue_annotation)

# Add tissue annotations to metadata
xenium_obj <- Seurat::AddMetaData(xenium_obj, metadata = manual_annotation)

# Remove cells that are annotated as 'remove'
xenium_obj <- subset(
    xenium_obj,
    subset = tolower(Tissue_annotation) != "remove"
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
