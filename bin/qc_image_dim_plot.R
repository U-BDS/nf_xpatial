#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package

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
        c("--dark_background"),
        type="character",
        default="F",
        help="Whether the plot has a dark background"),
    make_option(
        c("--group_by"),
        type="character",
        default=NULL,
        help="What to group the analysis by"),
    make_option(
        c("--width"),
        type="integer",
        default=2000,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=1500,
        help="Height of the plot"),
    make_option(
        c("--background"),
        type="character",
        default="transparent",
        help="Type of background for the plot"),
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
        c("--fov"),
        type="character",
        default=NULL,
        help="Name of FOV to plot"),
    make_option(
        c("--colors"),
        type="character",
        default="polychrome",
        help="Color palette for plot")
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

xenium.objs <- readRDS(
    file = opt$input
)

#################
#### DIM PLOT ###
#################

# TODO: How to do colors
# Need to get all the groups and assign them a color

# Check if the input was a list of objects or a single object
img_dim_plot <- NULL
if ( typeof(xenium.objs) != "list" ) {

    img_dim_plot <- 
        ImageDimPlot(
            xenium.objs,
            dark.background = opt$dark_background,
            group.by = opt$group_by,
            fov = opt$fov,
            cols = opt$colors
        ) + 
        ggtitle(xenium.objs@project.name) +
        coord_fixed() +
        scale_x_reverse()    

} else {

    fig_list = list()

    for (i in 1:length(xenium.objs)){
        fig_list[[i]] <- ImageDimPlot(
                xenium.objs[[i]],
                dark.background = opt$dark_background,
                group.by = opt$group_by,
                fov = opt$fov,
                cols = opt$colors
            ) + 
            ggtitle(xenium.objs[[i]]@project.name) +
            coord_fixed() +
            scale_x_reverse()
    }

    img_dim_plot <- wrap_plots(fig_list, nrow = opt$nrows, ncol = opt$ncols) + plot_layout(guides = "collect")
}

# Output the plot
png(
    opt$outfile,
    width = opt$width,
    height = opt$height,
    bg = opt$background
)

plot(img_dim_plot)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
