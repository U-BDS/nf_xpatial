#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(arrow)      # For reading/writing Parquet files
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(Seurat)     # Main analysis package

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
        c("-o", "--outfile"),
        type="character",
        default="outliers.csv",
        metavar="path",
        help="The csv file containing computed area quantiles"),
    make_option(
        c("--min_quantile"),
        type="double",
        default=0.25,
        help="The minimum quantile"),
    make_option(
        c("--max_quantile"),
        type="double",
        default=0.75,
        help="The maximum quantile")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide a Xenium object", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

##########################
### CALCULATE OUTLIERS ###
##########################

area_data <- as.numeric(xenium_obj@meta.data$Cell_Area)

# Calculate the quantiles

q1 <- quantile(area_data, opt$min_quantile)
q3 <- quantile(area_data, opt$max_quantile)

# interquantile range (iqr)
iqr <- q3 - q1

# boundaries
lower_bound <- q1 - (iqr * 1.5)
upper_bound <- q3 + (iqr * 1.5)

# Store outlier cutoffs in dataframe
outliers <- ifelse(area_data < lower_bound | area_data > upper_bound, TRUE, FALSE)
outlier_df <- data.frame(
    area = area_data,
    sample = xenium_obj$Sample,
    outlier = outliers,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    stringsAsFactors = FALSE
)

#################
### SAVE DATA ###
#################

# Save the outlier information to csv
write.csv(
    outlier_df,
    file = opt$outfile,
    row.names = TRUE,
    quote = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
