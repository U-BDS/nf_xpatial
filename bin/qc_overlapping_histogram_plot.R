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
        help="R Object to be analyzed"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="vln_plot.png",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("--width"),
        type="integer",
        default=1000,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=600,
        help="Height of the plot"),
    make_option(
        c("--nrows"),
        type="integer",
        default=NULL,
        help="Number of rows for the plot"),
    make_option(
        c("--ncols"),
        type="integer",
        default=1,
        help="Number of cols for the plot (if there are multiple samples)"),
    make_option(
        c("--alpha"),
        type="double",
        default=0.7,
        help="Alpha value for the points"),
    make_option(
        c("--fill"),
        type="character",
        default="skyblue",
        help="Fill color"),
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
        help="x-intercept for the plot")
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
    header = TRUE
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
png(
    opt$outfile,
    width = opt$width,
    height = opt$height
)

plot(overlapping_histogram_plot)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
