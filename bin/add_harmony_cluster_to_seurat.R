#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(Seurat)     # Main analysis package
library(stringr)    # For string manipulation

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
        c("-c", "--clusters"),
        type="character",
        default=NULL,
        metavar="path",
        help="The clusters csv to place on xenium object"),
    make_option(
        c("-e", "--embeddings"),
        type="character",
        default=NULL,
        metavar="path",
        help="The Embeddings csv to place on xenium object"),
    make_option(
        c("-l", "--loadings"),
        type="character",
        default=NULL,
        metavar="path",
        help="The Loadings csv to place on xenium object"),
    make_option(
        c("-s", "--stdev"),
        type="character",
        default=NULL,
        metavar="path",
        help="The Stdev csv to place on xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to keep during conversion"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="harmony_cluster_xenium_obj.rds",
        metavar="path",
        help="The output name for the seurat object")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Gather lists of files
cluster_file_list <- str_split_1(opt$clusters, pattern = ",")
embeddings_file_list <- str_split_1(opt$embeddings, pattern = ",")
loadings_file_list <- str_split_1(opt$loadings, pattern = ",")
stdev_file_list <- str_split_1(opt$stdev, pattern = ",")

# Set default assay
DefaultAssay(xenium_obj) <- opt$assay

##############################
### ADD REDUCED DIMENSIONS ###
##############################

print("Adding reduced dimensions...")

print("Parsing embeddings, loadings and stdev files")
print(embeddings_file_list)
print(loadings_file_list)
print(stdev_file_list)

# Parse file names for information
reduction_metadata_df <- t(data.frame(
    lapply(c(embeddings_file_list, loadings_file_list, stdev_file_list),
    function(x) {
        str_split_1(basename(x), pattern = "_")[c(1:4)]
    })
))

rownames(reduction_metadata_df) <- NULL
colnames(reduction_metadata_df) <- c("reduction", "dim", "res", "type")

print(reduction_metadata_df)

# Get the unique combinations of dimension, resolution, and reduction

reduction_metadata_df <- unique(reduction_metadata_df[, !(colnames(reduction_metadata_df) %in% c("type"))])
reduc_dim_res_list <- paste0(
    reduction_metadata_df[,"reduction"], 
    "_",
    reduction_metadata_df[,"dim"],
    "_",
    reduction_metadata_df[,"res"]
)
print(reduc_dim_res_list)

# Iterate through each unique combination, grab the corresponding files, and add to the Seurat object
for (reduc_dim_res in reduc_dim_res_list) {

    ######################
    ### ADD EMBEDDINGS ###
    ######################

    print(paste0("Grabbing files for  ",reduc_dim_res))

    # Obtain the embeddings file 
    print(paste0("Grabbing embeddings for  ",reduc_dim_res))
    embeddings_df <- read.csv(
        embeddings_file_list[grepl(reduc_dim_res, embeddings_file_list)][[1]],
        header = TRUE
    )

    rownames(embeddings_df) <- embeddings_df[,"Index"]
    embeddings_df <- embeddings_df[, !(colnames(embeddings_df) %in% c("Index"))]

    ####################
    ### ADD LOADINGS ###
    ####################

    # Should only be one file that matches
    loadings_file <- loadings_file_list[grepl(reduc_dim_res, loadings_file_list)]
    loadings_mtx <- matrix()
    if (length(loadings_file) > 1) {
        stop(paste0("More than one loadings file found for ", reduc_dim_res))

    } else if (length(loadings_file) == 1) {

        print(loadings_file[1])
        print(file.info(loadings_file[1])$size)
        if (file.info(loadings_file[1])$size <= 1) {
            print(paste0("Loadings file is empty for  ",reduc_dim_res))

        } else {
            print(paste0("Adding loadings for  ",reduc_dim_res))
            loadings_df <- read.csv(
                loadings_file[1],
                header = TRUE
            )
            rownames(loadings_df) <- loadings_df[,"Index"]
            loadings_df <- loadings_df[, !(colnames(loadings_df) %in% c("Index"))]

            loadings_mtx <-  as.matrix(loadings_df)
        }
    }

    #################
    ### ADD STDEV ###
    #################

    print(paste0("Grabbing stdev for  ",reduc_dim_res))
    stdev_df <- read.csv(
        stdev_file_list[grepl(reduc_dim_res, stdev_file_list)][[1]],
        header = TRUE
    )

    ############################
    ### CREATE DIM REDUCTION ###
    ############################
    print(paste0("Creating DimReducObject for  ",reduc_dim_res))
    xenium_obj[[reduc_dim_res]] <- CreateDimReducObject(
        embeddings = as.matrix(embeddings_df),
        loadings = loadings_mtx,
        stdev = as.numeric(stdev_df[[1]]),
        assay = opt$assay,
        key = paste0(reduc_dim_res, "_")
    )

}

############################
### ADD CLUSTER METADATA ###
############################

### TODO: see issue #21 for further evaluation on handling the cases where metadata may have missing cells

print("Adding cluster metadata...")
for (cluster_file in cluster_file_list) {

    if (file.exists(cluster_file)) {
        clusters_df <- read.csv(cluster_file, row.names = 1)
        
        # set clusters to be a sorted factor
        clusters_df[[1]] <- factor(clusters_df[[1]],
                                   levels = sort(unique(clusters_df[[1]])))

        # Ensure that the order of rownames matches between xenium metadata and new metadata
        common_cells <- intersect(rownames(xenium_obj@meta.data), rownames(clusters_df))
        
        # add warning in case common_cells are not matching all cells present in obj
        
        if (all(common_cells %in% rownames(xenium_obj@meta.data))) {
          print("Cell names between object and clustering metadata match. Adding cell cluster information")
          
          # Add the new metadata columns to the xenium object's metadata
          xenium_obj <- AddMetaData(object = xenium_obj,
                                    metadata = clusters_df)
          
        } else {
          warning("Cell names between object and clustering metadata DO NOT match. Skipping addition of cell cluster information!")
          
          # Subset both metadata to only include matching cells (see issue #21)
          #xenium_metadata <- xenium_obj@meta.data[common_cells, ]
          #clusters_df <- clusters_df[common_cells, ]
        }
      
    } else {
        print(paste0("File does not exist: ", cluster_file))
    }
}

#################
### SAVE DATA ###
#################

saveRDS(
    object = xenium_obj,
    file = paste0("d", opt$dim, "_r", opt$res, ".", opt$outfile)
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
