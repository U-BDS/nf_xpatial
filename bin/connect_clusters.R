#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(Seurat)     # Main analysis package
library(Banksy)
library(SummarizedExperiment) 

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
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to use"),
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
# Note this is a temporary object as input for connectClusters
xenium_obj_tmp <- readRDS(file = opt$input)

# Read in the banksy cluster info
metadata_df <- as.data.frame(read.table(
    file = opt$banksy_clust_info,
    header = TRUE,
    sep = ","
))

# Set the rownames and remove the index column
rownames(metadata_df) <- metadata_df$Index

metadata_df$Index <- NULL

##################################
### ADD PRE-CONNECTED METADATA ###
##################################

DefaultAssay(xenium_obj_tmp) <- opt$assay

## TODO: see issue #21 for further evaluation on handling cases where metadata has missing cells

# Ensure that the order of rownames matches between xenium metadata and new metadata
common_cells <- intersect(rownames(xenium_obj_tmp@meta.data), rownames(metadata_df))

## create a tmp Seurat object as input to 

if (all(common_cells %in% rownames(xenium_obj_tmp@meta.data))) {
  print("Cell names between object and clustering metadata match. Adding cell cluster information")
  
  # Add the new metadata columns to the xenium object's metadata
  xenium_obj_tmp <- AddMetaData(object = xenium_obj_tmp,
                                metadata = metadata_df)
  
  head(xenium_obj_tmp@meta.data)
  
} else {
  warning("Cell names between object and clustering metadata DO NOT match. Skipping addition of cell cluster information!")
  
  # Subset both metadata to only include matching cells (see issue #21)
  # xenium_metadata <- xenium_obj@meta.data[common_cells, ]
  # metadata_df <- metadata_df[common_cells, ]
  # metadata_df
}

###############################
### CONNECT BANKSY CLUSTERS ###
###############################

# Connect the clusters
xenium_obj_connected <- connectClusters(
    as.SingleCellExperiment(xenium_obj_tmp, assay = opt$assay)
)

# remove tmp data
rm(xenium_obj_tmp)
gc()

metadata_connected <- colData(xenium_obj_connected) %>% as.data.frame() %>% dplyr::select(!ident)

print("Head of post-connected metadata")
head(metadata_connected)

###################################
### ADD POST-CONNECTED METADATA ###
###################################

# re-read xenium object to add final metadata
### NOTE: this may create un-needed tmp memory ###
### If large, we should split and addition of final meta as separate processed ###
xenium_obj <- readRDS(file = opt$input)

DefaultAssay(xenium_obj) <- opt$assay

# Again, see issue #21 (for now if not equal, this step should crash)
# If filtering is needed, at we can leverage the cell names from xenium_obj_connected
stopifnot(
  all(rownames(xenium_obj@meta.data) == rownames(metadata_connected))
)

xenium_obj <- AddMetaData(object = xenium_obj,
                          metadata = metadata_connected)

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