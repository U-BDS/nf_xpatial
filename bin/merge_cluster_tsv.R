#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(arrow)      # For reading/writing Parquet files
library(dplyr)      # For data manipulation
library(jsonlite)   # For working with JSON data
library(knitr)      # For generating reports
library(optparse)   # For parsing commandline arguments
library(progressr)  # For progress bars
library(purrr)      # Functional programming tools
library(Seurat)     # Main analysis package

# Plotting
library(patchwork)  # For combining plots

# Set options
options(future.globals.maxSize = 16384 * 1024^2)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium results to be analyzed"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="merged_xenium_obj.rds",
        metavar="path",
        help="The output name for the xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default="Xenium",
        help="The assay name to be merged")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

###########################
### MERGE XENIUM OBJECT ###
###########################

merged_xenium_obj <- NULL

if (typeof(xenium_obj) == "list") {
    if (length(xenium_obj) >= 2) {
        merged_xenium_obj <- merge(
            x = xenium_obj[[1]],
            y = xenium_obj[2:length(xenium_obj)],
            merge.data = TRUE
        )
    } else {
        merged_xenium_obj <- xenium_obj[[1]]
    }

} else {
    merged_xenium_obj <- xenium_obj
}

DefaultAssay(merged_xenium_obj) <- opt$assay

if (opt$assay == "Xenium") {
    merged_xenium_obj <- JoinLayers(merged_xenium_obj)
}
#################
### SAVE DATA ###
#################

saveRDS(
    object = merged_xenium_obj,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
