#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package
library(stringr)

# Plotting
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(Polychrome)

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The input dataframe"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="overlap_hist_plot.png",
        metavar="path",
        help="The output name for the image"),
    make_option(
        c("--width"),
        type="integer",
        default=3000,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=3000,
        help="Height of the plot"),
    make_option(
        c("--alpha"),
        type="double",
        default=0.7,
        help="Alpha value for the points"),
    make_option(
        c("--color"),
        type="character",
        default="black",
        help="Color for the histogram"),
    make_option(
        c("--bins"),
        type="integer",
        default=30,
        help="Number of bins for the plot"),
    make_option(
        c("--threshold"),
        type="integer",
        default=1000,
        help="x-intercept (Area) for the plot")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide a dataframe as input.", call. = FALSE)
}

##################
### LOAD INPUT ###
##################

# Read in the xenium object
overlapping_histogram_df <- read.csv(
    file = opt$input,
    header = TRUE,
    row.names = 1
)

#######################
#### HISTOGRAM PLOT ###
#######################

overlapping_histogram_plot <- 
    ggplot(overlapping_histogram_df, aes(x = area, fill = sample)) +
    geom_histogram(
        position = "identity",
        alpha = opt$alpha,
        color = opt$color,
        bins = opt$bins
    ) +
    geom_vline(
        xintercept = opt$threshold,
        linetype = "dashed",
        color = "red",
        size = 1
    ) +
    geom_vline(
        aes(xintercept = upper_bound, color = sample),
        linetype = "dotted",
        size = 1) +
    labs(
        title = "Overlapping Histogram Plot",
        x = "Area",
        y = "Count"
    ) +
    theme_minimal() +
    theme(
        legend.position = "top"
    )

# Output the plot
ggsave(
    opt$outfile,
    plot = overlapping_histogram_plot,
    width = opt$width,
    height = opt$height,
    units = "px"
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
