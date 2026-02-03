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
        help="Fill color within binwidth of a histogram"),
    make_option(
        c("--color"),
        type="character",
        default="black",
        help="Color for the histogram"),
    make_option(
        c("--binwidth"),
        type="integer",
        default=1,
        help="binwidth of histogram")
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

# adjust ncols based on sample number if ncols <=1 (default)
# otherwise leave as user-selected number
if (opt$ncols <= 1) {
  if (length(xenium_objs) > 4) {
    opt$ncols <- ceiling(length(xenium_objs)/4) # round up beyond 4 samples
  }
}

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

###################
### OUTPUT PLOT ###
###################

# Calcuate width if not provided
indiv_plot_width <- 1500

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
    plot = histogram_plot,
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
