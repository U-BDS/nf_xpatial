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
        default="box_plot.png",
        metavar="path",
        help="The output name for the image"),
    make_option(
        c("--width"),
        type="integer",
        default=3500,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=2000,
        help="Height of the plot"),
    make_option(
        c("--alpha"),
        type="double",
        default=0.7,
        help="Alpha value for the boxplot")
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

# Read in the dataframe

box_plot_df <- read.csv(
    file = opt$input,
    header = TRUE,
    row.names = 1
)

#################
#### BOX PLOT ###
#################

box_plot <- 
    ggplot(
        data = box_plot_df, 
        aes(x = sample, y = area, fill = sample)
    ) +
    geom_boxplot(
        alpha = opt$alpha,
        outlier.shape = NA,
    ) +
    geom_point(
        data = subset(box_plot_df, outlier == TRUE),
        aes(x = sample, y = area),
        color = "red",
        size = 2,
        shape = 16
    ) +
    labs(title = "Box Plot",
            x = "Sample",
            y = "Area") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

# Output the plot
ggsave(
    opt$outfile,
    plot = box_plot,
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
