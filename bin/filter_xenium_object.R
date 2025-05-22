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
library(stringr)   # For string manipulation

# Set options
#options(future.globals.maxSize = 8192 * 1024^2)

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
        c("-o", "--outfile"),
        type="character",
        default="filtered_xenium_obj.rds",
        metavar="path",
        help="The filtered xenium object"),
    make_option(
        c("--skip_col_removal"),
        action="store_true",
        default=FALSE,
        help="Skip the removal of columns from the metadata"),
    make_option(
        c("-c", "--columns"),
        type="character",
        default="nCount_BlankCodeword,nFeature_BlankCodeword,nCount_ControlCodeword,nFeature_ControlCodeword,nCount_ControlProbe,nFeature_ControlProbe",
        help="The list of columns to be removed from the metadata"),
    make_option(
        c("--skip_percentile_filter"),
        action="store_true",
        default=FALSE,
        help="Skip filtering data by percentile"),
    make_option(
        c("--min_percentile"),
        type="integer",
        default=-1,
        help="The minimum percentile to filter by"),
    make_option(
        c("--max_percentile"),
        type="integer",
        default=101,
        help="The maximum percentile to filter by"),
    make_option(
        c("--skip_nFeature_filter"),
        action="store_true",
        default=FALSE,
        help="Skip filtering data by nFeature"),
    make_option(
        c("--min_nFeature"),
        type="integer",
        default=NULL,
        help="The minimum number of nFeatures to filter by"),
    make_option(
        c("--max_nFeature"),
        type="integer",
        default=NULL,
        help="The maximum number of nFeatures to filter by"),
    make_option(
        c("--skip_nCount_filter"),
        action="store_true",
        default=FALSE,
        help="Skip filtering data by nFeature"),
    make_option(
        c("--min_nCount"),
        type="integer",
        default=NULL,
        help="The minimum number of nCounts to filter by"),
    make_option(
        c("--max_nCount"),
        type="integer",
        default=NULL,
        help="The maximum number of nCounts to filter by")
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

###################
### FILTER DATA ###
###################

# Initialize variables to be output to assess filtering
pre_filtered_cells <- sum(table(xenium_obj@meta.data$orig.ident))
percentile_filtered <- 0
nFeature_filtered <- 0
nCount_filtered <- 0

if (!opt$skip_col_removal){
    # Remove unnecessary columns from metadata
    cols_to_remove <- str_split_1(opt$columns, ",")

    for (col in cols_to_remove) {
        if (col %in% colnames(xenium_obj@meta.data)){
            xenium_obj@meta.data[[col]] <- NULL
        }
    }
}

# Remove control assays

# Remove features

# filter by percentile
if (!opt$skip_percentile_filter){
    # Save prefiltering count
    pre_percentile_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))

    # Calculate percentiles
    percentiles <- data.frame(percentile = seq(from = 0, to = 100, by = 10))
    percentiles$count <- matrixStats::colQuantiles(
        as.matrix(xenium_obj$nFeature_Xenium),
        probs = percentiles$percentile / 100
    )

    # Get cutoff values
    min_count <- percentiles[percentiles$percentile == opt$min_percentile, "count"]
    max_count <- percentiles[percentiles$percentile == opt$max_percentile, "count"]

    # Filter my minimum counts
    tryCatch(
        {
            xenium_obj <- subset(xenium_obj, nFeature_Xenium >= min_count)
        }, error = function(e) {
            message("Error in min percentile filtering: ", e)
        }
    )

    # Filter by maximum counts
    tryCatch(
        {
            xenium_obj <- subset(xenium_obj, nFeature_Xenium <= max_count)
        }, error = function(e) {
            message("Error in max percentile filtering: ", e)
        }
    )

    # Calculate cells removed
    post_percentile_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))
    percentile_filtered <- pre_percentile_filter_count - post_percentile_filter_count
}

# filter by nCount
if (!opt$skip_nCount_filter){
    # Save prefiltering count
    pre_nCount_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))

    # Filter by minimum nFeature
    if (!is.null(opt$min_nCount)){
        xenium_obj <- subset(xenium_obj, nCount_Xenium >= opt$min_nCount)
    }

    # Filter by maximum nFeature
    if (!is.null(opt$max_nCount)){
        xenium_obj <- subset(xenium_obj, nCount_Xenium <= opt$max_nCount)
    }

    # Calculate cells removed
    post_nCount_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))
    nCount_filtered <- pre_nCount_filter_count - post_nCount_filter_count
}

# filter by nFeature
if (!opt$skip_nFeature_filter){
    # Save prefiltering count
    pre_nFeature_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))

    # Filter by minimum nFeature
    if (!is.null(opt$min_nFeature)){
        xenium_obj <- subset(xenium_obj, nFeature_Xenium >= opt$min_nFeature)
    }

    # Filter by maximum nFeature
    if (!is.null(opt$max_nFeature)){
        xenium_obj <- subset(xenium_obj, nFeature_Xenium <= opt$max_nFeature)
    }

    # Calculate cells removed
    post_nFeature_filter_count <- sum(table(xenium_obj@meta.data$orig.ident))
    nFeature_filtered <- pre_nFeature_filter_count - post_nFeature_filter_count
}

#################
### SAVE DATA ###
#################

# Save the filtered xenium object
saveRDS(
    object = xenium_obj,
    file = opt$outfile 
)

# Save filtering stats as csv
post_filtered_cells <- sum(table(xenium_obj@meta.data$orig.ident))
cell_numbers <- data.frame(
    sample_name = xenium_obj@project.name,
    cells_pre_filter = pre_filtered_cells,
    cells_post_filter = post_filtered_cells,
    cells_filtered_percentile = percentile_filtered,
    cells_filtered_nFeature = nFeature_filtered,
    cells_filtered_nCount = nCount_filtered
)

write.csv(
    cell_numbers,
    file = paste0(
        gsub(".rds", "", opt$outfile),
        "_filtering_stats.csv"
    ),
    row.names = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
