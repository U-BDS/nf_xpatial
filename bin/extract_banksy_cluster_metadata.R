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

################################
### EXTRACT CLUSTER METADATA ###
################################

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

col_data <- colData(spe_xenium_obj)

clusts <- col_data[ grepl("^clust_BSKY", colnames(col_data)) ]
clusts$Index <- rownames(clusts)
clusts <- clusts[rev(colnames(clusts))]

###################
### WRITE TABLE ###
###################

# Save the cluster_table
write.table(
    clusts,
    file = opt$outfile,
    sep = "\t",
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
