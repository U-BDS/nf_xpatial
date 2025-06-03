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
        c("--binwidth"),
        type="integer",
        default=1,
        help="Number of cols for the plot (if there are multiple samples)")
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

#######################
#### HISTOGRAM PLOT ###
#######################

# Check if the input was a list of objects or a single object
histogram_plot <- NULL
if ( typeof(xenium_objs) != "list" ) {
    histogram_plot <- 
        ggplot(
            data = data.frame(Cell_Area = xenium_objs@meta.data$Cell_Area), 
            aes(x = Cell_Area)
        ) +
        geom_histogram(binwidth = opt$binwidth, fill = opt$fill, color = opt$color, alpha = opt$alpha) +
        labs(
            title = paste0("Histogram of Cell Area for ", xenium_objs$Sample),
            x = "Area",
            y = "Frequency") +
        theme_minimal()

} else {

    fig_list = list()

    for (i in 1:length(xenium_objs)){
        fig_list[[i]] <- 
            ggplot(
                data = data.frame(Cell_Area = xenium_objs[[i]]@meta.data$Cell_Area), 
                aes(x = Cell_Area)
            ) +
            geom_histogram(binwidth = opt$binwidth, fill = opt$fill, color = opt$color, alpha = opt$alpha) +
            labs(
                title = paste0("Histogram of Cell Area for ", xenium_objs[[i]]$Sample),
                x = "Area",
                y = "Frequency") +
            theme_minimal()
    }

    histogram_plot <- wrap_plots(fig_list, nrow = opt$nrows, ncol = opt$ncols)
}

# Output the plot
png(
    opt$outfile,
    width = opt$width,
    height = opt$height
)

plot(histogram_plot)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
