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
        type="character",
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

cluster_metadata <- NULL

if (class(xenium_obj) == "Seurat") {
    # Set default assay
    DefaultAssay(xenium_obj) <- opt$assay

    # Extract cluster data
    cluster_metadata <- xenium_obj@meta.data[grepl("seurat_clusters", colnames(xenium_obj@meta.data))]

    # Rename seurat cluster column with a similar prefix to BANKSY
    colnames(cluster_metadata) <- gsub("seurat_clusters", "clust_HMY_", colnames(cluster_metadata))

    # Add dim and res to column name
    colnames(cluster_metadata) <- paste0(colnames(cluster_metadata), opt$param_string)

    # Create a column for the cell ids
    cluster_metadata$Index <- rownames(cluster_metadata)
    rownames(cluster_metadata) <- NULL

    # Rearrange columns
    cluster_metadata <- cluster_metadata[, c(2,1)]

} else if (class(xenium_obj) == "SpatialExperiment") {
    # Set default assay
    mainExpName(xenium_obj) <- opt$assay
    col_data <- colData(xenium_obj)

    cluster_metadata <- col_data[ grepl("^clust_BSKY", colnames(col_data)) ]
    cluster_metadata$Index <- rownames(cluster_metadata)
    cluster_metadata <- cluster_metadata[rev(colnames(cluster_metadata))]

} else {
    stop("Input object is not of class SeuratObject or SpatialExperiment.")
}

###################
### WRITE TABLE ###
###################

# Save the cluster_table
write.table(
    cluster_metadata,
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
