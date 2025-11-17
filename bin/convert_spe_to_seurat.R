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
library(Banksy)
library(SummarizedExperiment) 

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-x", "--xenium"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium object to add spe information onto"),
    make_option(
        c("-s", "--spe_obj"),
        type="character",
        default=NULL,
        metavar="path",
        help="The spe object to extract information from"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to keep during conversion"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="spe_to_xenium_obj.rds",
        metavar="path",
        help="The filtered xenium object")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$xenium)
# Set the default assay
DefaultAssay(xenium_obj) <- opt$assay

# Read in spe object
spe_obj <- readRDS(file = opt$spe_obj)
# Set the default assay
mainExpName(spe_obj) <- opt$assay

################################
### EXTRACT CLUSTER METADATA ###
################################

col_data <- colData(spe_obj)

clusts <- col_data[ grepl("^clust_BSKY", colnames(col_data)) ]
colnames( clusts) <- c("seurat_clusters")
clusts <- clusts[rev(colnames(clusts))]

############################
### ADD CLUSTER METADATA ###
############################

## TODO: see issue #21 for further evaluation on handling cases where metadata has missing cells
# Ensure that the order of rownames matches between xenium metadata and new metadata
common_cells <- intersect(rownames(xenium_obj@meta.data), rownames(clusts))

## create a tmp Seurat object as input to 

if (all(common_cells %in% rownames(xenium_obj@meta.data))) {
    print("Cell names between object and clustering metadata match. Adding cell cluster information")

    # Add the new metadata columns to the xenium object's metadata
    xenium_obj <- AddMetaData(
        object = xenium_obj,
        metadata = as.data.frame(clusts)
    )

    head(xenium_obj@meta.data)

} else {
    warning("Cell names between object and clustering metadata DO NOT match. Skipping addition of cell cluster information!")
}

############################
### EXTRACT REDUCED DIMS ###
############################

# Save each sample-specific data frame to a csv file
for (dim_name in reducedDimNames(spe_obj)) {
    reduced_dim_data <- reducedDim(spe_obj, dim_name)

    converted_dim_name <- dim_name
    if (dim_name == "UMAP_BANKSY_harmony") {
        converted_dim_name <- "umap"

        colnames(reduced_dim_data) <- c("umap_1", "umap_2")
    } else {
        converted_dim_name <- gsub("BANKSY_", "", dim_name)
    }

    print(reduced_dim_data)

    xenium_obj[[converted_dim_name]] <- CreateDimReducObject(
        embeddings = as.matrix(reduced_dim_data),
        assay = opt$assay,
        key = paste0(converted_dim_name, "_")
    )
}

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
