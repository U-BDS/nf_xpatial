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
        default=NULL,
        help="The assay to operate on"),
    make_option(
        c("-r", "--res"),
        type="integer",
        help="The resolution to use for Find clusters"
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="find_clusters_xenium_obj.rds",
        metavar="path",
        help="The output name for the xenium object"))

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}

if (is.null(opt$assay)) {
    print_help(opt_parser)
    stop("Please provide the assay to integrate on as input.", call. = FALSE)
}

if (is.null(opt$reduction)) {
    print_help(opt_parser)
    stop("Please provide the reduction name to integrate on as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

######################
### FIND NEIGHBORS ###
######################

DefaultAssay(xenium_obj) <- opt$assay
xenium_obj <- FindClusters (
    xenium_obj,
    dims = 1:opt$resolution
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
