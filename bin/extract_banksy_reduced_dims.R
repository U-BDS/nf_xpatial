#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)          # For parsing commandline arguments
library(Seurat)            # For handling Seurat objects
library(SpatialExperiment) # For handling SpatialExperiment objects
library(stringr)           # For string manipulation

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
        help="The assay to keep during conversion"),
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

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

if (length(reducedDimNames(spe_xenium_obj)) == 0) {
    stop("No reduced dimensions found in the Spatial Experiment object")
}

############################
### EXTRACT REDUCED DIMS ###
############################

# Save each sample-specific data frame to a csv file
for (dim_name in reducedDimNames(spe_xenium_obj)) {
    reduced_dim_data <- reducedDim(spe_xenium_obj, dim_name)

    reduced_dim_df <- data.frame(
        Index = rownames(reduced_dim_data),
        reduced_dim_data,
        stringsAsFactors = FALSE
    )

    write.csv(
        reduced_dim_df,
        file = paste0(dim_name, "_", opt$outfile),
        row.names = FALSE
    )
}

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
