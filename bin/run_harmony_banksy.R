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
        help="The xenium object to be filtered"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        metavar="path",
        help="The assay to keep during conversion"),
    make_option(
        c("--max_iter"),
        type="int",
        default=60,
        help="Maximum number of iterations for Banksy clustering"),
    make_option(
        c("--ncores"),
        type="int",
        default=1,
        help="Number of cores to use for Harmony"),
   make_option(
        c("--nPCs"),
        type="integer",
        default=NULL,
        help="nPCs values for Banksy clustering"),
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

##############
### BANKSY ###
##############

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

# Compute banksy matrix
harmony_embedding <- RunHarmony(
    data_mat = reducedDim(spe_xenium_obj, "BANKSY_pca")[, 1:opt$nPCs],
    meta_data = colData(spe_xenium_obj),
    vars_use = c("orig.ident"),
    plot_convergence = FALSE,
    max_iter = opt$max_iter,
    verbose = TRUE,
    ncores = opt$ncores
)

reducedDim(spe_xenium_obj, "BANKSY_harmony") <- harmony_embedding

#################
### SAVE DATA ###
#################

# Save the filtered xenium object
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
