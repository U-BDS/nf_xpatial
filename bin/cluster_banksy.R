#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(Banksy)           # For Banksy clustering
library(harmony)          # For data integration
library(optparse)          # For parsing commandline arguments
library(Seurat)            # For handling Seurat objects
library(SpatialExperiment) # For handling SpatialExperiment objects
library(tidyverse)         # For data manipulation

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
        c("--use_harmony"),
        action="store_true",
        default=FALSE,
        help="Whether to use harmony for clustering"),
    make_option(
        c("--param_str"),
        type="character",
        default=NULL,
        metavar="path",
        help="The string of parameters to use"),
    make_option(
        c("--lambda"),
        type="double",
        default=NULL,
        help="lambda values for Banksy clustering"),
   make_option(
        c("--nPCs"),
        type="integer",
        default=NULL,
        help="nPCs values for Banksy clustering"),
    make_option(
        c("--res"),
        type="double",
        default=NULL,
        help="Resolution values for Banksy clustering"),
    make_option(
        c("--agf"),
        action="store_true",
        default=FALSE,
        help="Whether to use the adaptive Gaussian filter (AGF) in Banksy"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="filtered_xenium_obj.rds",
        metavar="path",
        help="The filtered xenium object")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
spe_xenium_obj <- readRDS(file = opt$input)

######################
### CLUSTER_BANKSY ###
######################

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

# Perform Banksy clustering
if (opt$use_harmony) {
    spe_xenium_obj <- clusterBanksy(
        spe_xenium_obj,
        dimred = "BANKSY_harmony",
        use_agf = opt$agf,
        resolution = opt$res,
        lambda = opt$lambda,
        ndims = opt$nPCs,
        seed = 1234,
        verbose = TRUE
  )
} else {
    spe_xenium_obj <- clusterBanksy(
        spe_xenium_obj,
        dimred = "BANKSY_pca",
        use_agf = opt$agf,
        resolution = opt$res,
        lambda = opt$lambda,
        ndims = opt$nPCs,
        seed = 1234,
        verbose = TRUE
  )
}

# add new column name with Banksy parameters containing cluster assignments
colData(spe_xenium_obj)[[paste0("clust_", opt$param_str)]] <- colData(spe_xenium_obj)[[ncol(colData(spe_xenium_obj))]]

#################
### SAVE DATA ###
#################

saveRDS(
    object = spe_xenium_obj,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
