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
xenium_obj <- readRDS(file = opt$input)

######################
### CONVERT TO SPE ###
######################

# Set the default assay
DefaultAssay(xenium_obj) <- opt$assay

# counts <- Seurat::GetAssayData(
#     xenium_obj,
#     slot = "counts"
# )

# metadata <- xenium_obj@meta.data

# spatial_coords <- data.frame (x = xenium_obj$x, y = xenium_obj$y)
# colData <- S4Vectors::DataFrame(metadata)

# spe <- SpatialExperiment::SpatialExperiment(
#     assays = list(counts = counts),
#     colData = colData
# )

# head(xenium_obj@meta.data)
# dim(xenium_obj)
# dim(spatial_coords)
# dim(spe)

# spatial_exp <- SpatialExperiment::SpatialExperiment(
#     assays = list(counts = counts),
#     colData = colData,
#     spatialCoords = as.matrix(spatial_coords)
# )

# Extract counts
single_cell_exp <- Seurat::as.SingleCellExperiment(
    xenium_obj,
    assay = opt$assay
)
spatial_coords <- data.frame (x = xenium_obj$x, y = xenium_obj$y) 

assays(single_cell_exp)

# Convert to spatial experiment object
spatial_exp <- SpatialExperiment(
    assays = assays(single_cell_exp),
    rowData = rowData(single_cell_exp),
    colData = colData(single_cell_exp),
    metadata = metadata(single_cell_exp),
    reducedDims = reducedDims(single_cell_exp),
    altExps = altExps(single_cell_exp),
    spatialCoords = as.matrix(spatial_coords)
)

# Set the default assay on the spatial experiment object
mainExpName(spatial_exp) <- opt$assay

#################
### SAVE DATA ###
#################

# Save the filtered xenium object
saveRDS(
    object = spatial_exp,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
