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
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to use for plotting (if applicable)"),
    make_option(
        c("--width"),
        type="integer",
        default=0,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=0,
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
    if (!is.null(opt$assay)){
        DefaultAssay(xenium_objs) <- opt$assay
    }

    img_feature_plot <- 
        ImageFeaturePlot(
            xenium_objs,
            features = features,
            cols = c("#000000","#FF0000")
        ) + 
        ggtitle(xenium_objs@project.name) +
        coord_fixed() +
        scale_x_reverse()    

} else {

    fig_list = list()

    for (i in 1:length(xenium_objs)){
        if (!is.null(opt$assay)){
            DefaultAssay(xenium_objs[[i]]) <- opt$assay
        }

        fig_list[[i]] <- ImageFeaturePlot(
            xenium_objs[[i]],
            features = features,
            cols = c("#000000","#FF0000")
        ) + 
        ggtitle(xenium_objs[[i]]@project.name) +
        coord_fixed() +
        scale_x_reverse()   
    }

    img_feature_plot <- wrap_plots(fig_list, nrow = opt$nrows, ncol = opt$ncols)
}

###################
### OUTPUT PLOT ###
###################

# Calcuate width if not provided
indiv_plot_width <- 1000

total_plot_width <- opt$width
if ( opt$width <= 0 ) {
    if (typeof(xenium_objs) != "list") {
        total_plot_width <- indiv_plot_width
    } else {
        col_number <- ifelse(opt$ncols > length(xenium_objs), length(xenium_objs), opt$ncols)
        total_plot_width <- indiv_plot_width * col_number
    }
}

# Calculate height if not provided
indiv_plot_height <- 1000

total_plot_height <- opt$height
if ( opt$height <= 0 ) {
    if (typeof(xenium_objs) != "list") {
        total_plot_height <- indiv_plot_height
    } else {
        total_plot_height <- indiv_plot_height * ceiling(length(xenium_objs) / opt$ncols)
    }
}

# Output the plot
ggsave(
    opt$outfile,
    plot = img_feature_plot,
    width = total_plot_width,
    height = total_plot_height,
    units = "px",
    limitsize = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
