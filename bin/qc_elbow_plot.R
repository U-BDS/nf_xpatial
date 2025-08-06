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
        c("-a", "--assay"),
        type="character",
        default="Xenium",
        help="The assay name"),
    make_option(
        c("-r", "--reduction_name"),
        type="character",
        default="pca",
        help="The reduction name to evaluate on the xenium object"),
    make_option(
        c("-n", "--ndims"),
        type="integer",
        default=50,
        help="The number of dims to generate"),
    make_option(
        c("--width"),
        type="integer",
        default=800,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=480,
        help="Height of the plot")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide a xenium object as input.", call. = FALSE)
}

##################
### LOAD INPUT ###
##################

# Read in the xenium object
xenium_objs <- readRDS(
    file = opt$input
)

###################
#### ELBOW PLOT ###
###################

elbow_plot <- 
    ElbowPlot(
        object = xenium_objs,
        ndims = opt$ndims,
        reduction = opt$reduction_name
    )

# Output the plot
png(
    opt$outfile,
    width = opt$width,
    height = opt$height
)

plot(elbow_plot)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
