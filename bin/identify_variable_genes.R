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
        help="The xenium results to be analyzed"),
    make_option(
        c("-p", "--percent"),
        type="integer",
        default=10,
        help="The percent of highly variable genes to select"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="hvg_list.csv",
        metavar="path",
        help="The output name for the gene list"))

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

#################
### NORMALIZE ###
#################

num_hvg <- as.integer(opt$percent / 100 * dim(xenium_obj)[1])

xenium_obj <- FindVariableFeatures(
    xenium_obj,
    select.method = "vst",
    nfeatures = num_hvg
)

hvg_list <- VariableFeatures(xenium_obj)

#################
### SAVE DATA ###
#################

write.table(
    hvg_list,
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
