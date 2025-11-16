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
        help="The xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to use"),
    make_option(
        c("-p", "--param_string"),
        type="string",
        help="The parameter string to use to identify clustering run"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="cluster_metadata.csv",
        metavar="path",
        help="The csv containing seurat cluster metadata")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

################################
### EXTRACT CLUSTER METADATA ###
################################

# Set default assay
DefaultAssay(xenium_obj) <- opt$assay

# Extract cluster data
clusts <- xenium_obj@meta.data[grepl("seurat_clusters", colnames(xenium_obj@meta.data))]

# Add dim and res to column name
colnames(clusts) <- paste0(
    colnames(clusts),
    opt$param_string
)

# Create a column for the cell ids
clusts$Index <- rownames(clusts)
rownames(clusts) <- NULL

# Rearrange columns
clusts <- clusts[, c(2,1)]

###################
### WRITE TABLE ###
###################

# Save the cluster_table
write.table(
    clusts,
    file = opt$outfile,
    sep = ",",
    quote = FALSE,
    row.names = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
