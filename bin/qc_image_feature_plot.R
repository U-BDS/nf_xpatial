#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package
library(stringr)    # String manipulation

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
        default="image_dim_plot.png",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("--width"),
        type="integer",
        default=500,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=1000,
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
        help="Number of cols for the plot"),
    make_option(
        c("--features"),
        type="character",
        default=NULL,
        help="Features to plot")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide a seurat object as input.", call. = FALSE)
}

##################
### LOAD INPUT ###
##################

# Read in the xenium object
xenium_objs <- readRDS(
    file = opt$input
)

# Parse features
features <- str_split_1(opt$features, ",")

###########################
#### IMAGE FEATURE PLOT ###
###########################

# TODO: How to do colors
# Need to get all the groups and assign them a color

# Check if the input was a list of objects or a single object
img_feature_plot <- NULL
if ( typeof(xenium_objs) != "list" ) {

    img_feature_plot <- 
        ImageFeaturePlot(
            xenium_objs,
            features = features
        ) + 
        ggtitle(xenium_objs@project.name) +
        coord_fixed() +
        scale_x_reverse()    

} else {

    fig_list = list()

    for (i in 1:length(xenium_objs)){
        fig_list[[i]] <- ImageFeaturePlot(
            xenium_objs[[i]],
            features = features
        ) + 
        ggtitle(xenium_objs[[i]]@project.name) +
        coord_fixed() +
        scale_x_reverse()   
    }

    img_feature_plot <- wrap_plots(fig_list, nrow = opt$nrows, ncol = opt$ncols)
}

# Output the plot
png(
    opt$outfile,
    width = opt$width,
    height = opt$height
)

plot(img_feature_plot)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
