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
        c("-o", "--outfile"),
        type="character",
        default="extracted_metadata.csv",
        metavar="path",
        help="The extracted metadata")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
spe_xenium_obj <- readRDS(file = opt$input)

###################################################
### EXTRACT BANKSY METADATA for XENIUM EXPLORER ###
###################################################

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

# Extract colData from the SpatialExperiment object
col_data <- as.data.frame(colData(spe_xenium_obj))

# Check if any column starts with clust_BSKY
clust_columns <- grep("^clust_BSKY", colnames(col_data), value = TRUE)
if (length(clust_columns) == 0) {
    stop("No columns starting with 'clust_BSKY' found in colData.")
}

# Check if 'Sample' and 'Cell_ID' columns exist
if (!"Sample" %in% colnames(col_data)) {
    stop("'Sample' column not found in colData.")
}

if (!"Cell_ID" %in% colnames(col_data)) {
    stop("'Cell_ID' column not found in colData.")
}

# Extract sample and cell_id columns
sample_column <- col_data$Sample
cell_id_column <- col_data$Cell_ID

# Loop through each clustering column and export data frames by sample
for (clust_column in clust_columns) {
    # Renamed the current clustering column to 'group'
    metadata_df <- data.frame(
        sample = sample_column,
        cell_id = cell_id_column,
        group = col_data[[clust_column]],
        stringsAsFactors = FALSE
    )
}

# split data by sample
grouped_data <- split(metadata_df[, c("cell_id", "group")], metadata_df$sample)

###################
### WRITE TABLE ###
###################

# Save each sample-specific data frame to a csv file
for (sample_name in names(grouped_data)) {
    write.csv(
        grouped_data[[sample_name]],
        file = paste0(sample_name, "_", opt$outfile),
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
