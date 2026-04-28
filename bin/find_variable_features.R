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
        c("-n", "--nfeatures"),
        type="integer",
        default=2000,
        help="The number of variable features to select"
    ),
    make_option(
        c("-p", "--percent"),
        type="double",
        default=0,
        help="The percent of highly variable genes to select"
    ),
    make_option(
        c("-s","--selection_method"),
        type="character",
        default="vst",
        help="The method to use for variable feature selection."
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="vf_xenium",
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

if (opt$percent > 0 && opt$nfeatures > 0) {
    stop("Please only specify -p/--percent or -n/--nfeatures, not both.", call. = FALSE)

} else if (opt$percent == 0 && opt$nfeatures == 0) {
    stop("Please specify either -p/--percent or -n/--nfeatures.", call. = FALSE)

}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

DefaultAssay(xenium_obj) <- opt$assay

##############################
### FIND VARIABLE FEATURES ###
##############################

num_features <- 0
if (opt$percent > 0) {
    num_features <- as.integer(opt$percent / 100 * dim(xenium_obj)[1])

} else if (opt$nfeatures > 0) {
    num_features <- opt$nfeatures

}

if (num_features > nrow(xenium_obj[[opt$assay]]$counts)) {
    warning(paste0(
        "The number of total features available in the current assay (",
        nrow(xenium_obj[[opt$assay]]$counts),
        ") is less than ", 
        opt$nfeatures,
        "\nFindVariableFeatures will result in using all available features."
        )
    )

    num_features <- nrow(xenium_obj[[opt$assay]]$counts)
}

xenium_obj <- FindVariableFeatures(
    object = xenium_obj,
    assay = opt$assay,
    selection.method = opt$selection_method,
    nfeatures = num_features
)

vf_list <- VariableFeatures(xenium_obj)

#################
### SAVE DATA ###
#################

saveRDS(
    object = xenium_obj,
    file = paste0(opt$outfile,".rds")
)

write.table(
    vf_list,
    file = paste0(opt$outfile,".csv"),
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
