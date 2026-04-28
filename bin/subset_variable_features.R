#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(arrow)             # For reading/writing Parquet files
library(dplyr)             # For data manipulation
library(optparse)          # For parsing commandline arguments
library(Seurat)            # Main analysis package
library(SpatialExperiment) # For handling SpatialExperiment objects

# Plotting
library(patchwork)  # For combining plots

# Set options
options(future.globals.maxSize = 16384 * 1024^2)

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
        c("-o", "--outfile"),
        type="character",
        default="xenium_obj_vf_subset.rds",
        metavar="path",
        help="The output name for the xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default="Xenium",
        help="The assay name to be merged"),
    make_option(
        c("-w", "--whitelist"),
        type="character",
        default="",
        metavar="path",
        help="The list of features/genes to keep in addition to variable features")
    )

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

# Set the default assay on the xenium object
DefaultAssay(xenium_obj) <- opt$assay

################################
### SUBSET VARIABLE FEATURES ###
################################

features_keep <- sort(VariableFeatures(xenium_obj))

# Read in the whitelist of features to keep, if provided
if (opt$whitelist != "") {
    features_whitelist <- readLines(opt$whitelist)

    features_keep <- sort(unique(c(features_keep, features_whitelist)))
}

xenium_obj <- subset(xenium_obj, features = features_keep)

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
