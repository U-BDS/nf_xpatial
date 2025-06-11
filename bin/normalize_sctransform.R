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
        c("-v", "--vst_flavor"),
        type="character",
        default="v2",
        help="The vst flavor to use"),
    make_option(
        c("-m", "--method"),
        type="character",
        default="glmGamPoi",
        help="The method to use"),
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

#################
### NORMALIZE ###
#################

xenium_obj <- SCTransform(
    xenium_obj,
    assay = opt$assay,
    return.only.var.genes = FALSE,
    verbose = TRUE,
    vst.flavor = opt$vst_flavor
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
