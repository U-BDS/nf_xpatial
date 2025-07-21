#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(Seurat)     # Main analysis package

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
        help="The xenium object"),
    make_option(
        c("-b", "--banksy_clust_info"),
        type="character",
        default=NULL,
        metavar="path",
        help="The BANKSY cluster information file"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="xenium_obj.rds",
        metavar="path",
        help="The output name for the seurat object")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium object as input.", call. = FALSE)
}

if (is.null(opt$banksy_clust_info)) {
    print_help(opt_parser)
    stop("Please provide the banksy cluster file to add to the xenium object.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Read in the banksy cluster info
metadata_df <- as.data.frame(read.table(
    file = opt$banksy_clust_info,
    header = TRUE,
    sep = ","
))

# Set the rownames and remove the index column
rownames(metadata_df) <- metadata_df$Index

metadata_df$Index <- NULL

####################
### ADD METADATA ###
####################

# Ensure that the order of rownames matches between xenium metadata and new metadata
common_cells <- intersect(rownames(xenium_obj@meta.data), rownames(metadata_df))

# Subset both metadata to only include matching cells
xenium_metadata <- xenium_obj@meta.data[common_cells, ]
metadata_df <- metadata_df[common_cells, ]

# Add the new metadata columns to the xenium object's metadata
xenium_obj@meta.data <- cbind(xenium_metadata, metadata_df)


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
