#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(Banksy)            # For Banksy clustering
library(optparse)          # For parsing commandline arguments
library(Seurat)            # For handling Seurat objects
library(SeuratWrappers)    # For Seurat wrappers around other tools

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to use"),
    make_option(
        c("--lambda"),
        type="double",
        default=NULL,
        help="lambda values for Banksy clustering"),
   make_option(
        c("--k_geom"),
        type="integer",
        default=NULL,
        help="Number of neighbors to use for Banksy"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="filtered_xenium_obj.rds",
        metavar="path",
        help="The filtered xenium object"),
    make_option(
        c("-f", "--features"),
        type="character",
        default="variable",
        help="The features to use for clustering (default: all, options: all, variable)")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Set the default assay on the xenium object
DefaultAssay(xenium_obj) <- opt$assay

######################
### RUN_BANKSY ###
######################

xenium_obj <- RunBanksy(
    xenium_obj,
    lambda = opt$lambda,
    k_geom = opt$k_geom,
    assay = opt$assay,
    dimx = "x",
    dimy = "y",
    slot = "data",
    features = opt$features,
    group = "Sample",
    verbose = TRUE,
    split.scale = FALSE,
    assay_name = paste0(opt$assay, "_BANKSY")
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
