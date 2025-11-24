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
library(pals)

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
        c("-g", "--gene_list"),
        type="character",
        default=NULL,
        metavar="path",
        help="The marker gene list"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="harmony_plots.png",
        metavar="path",
        help="The output name for the harmony_plot"),
    make_option(
        c("--width"),
        type="integer",
        default=2500,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=3000,
        help="Height of the plot"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default="Xenium",
        help="The assay name to be merged"),
    make_option(
        c("-c", "--cluster_col"),
        type="character",
        default="seurat_clusters",
        help="The column name to pull for cluster numbers")
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
xenium_obj <- readRDS(
    file = opt$input
)
DefaultAssay(xenium_obj) <- opt$assay

# Read in marker gene list
marker_gene_list <- read.csv(opt$gene_list, sep = ",", header=FALSE)
marker_gene_list <- unlist(marker_gene_list, use.names = FALSE)

# Set cluster column
cluster_col <- opt$cluster_col
print(cluster_col)

#################
#### VLN PLOT ###
#################

# Get colors
cols <- NULL
if (length(unique(xenium_obj@meta.data[[cluster_col]])) <= 32) {
    cols <- as.vector(pals::glasbey(length(unique(xenium_obj@meta.data[[cluster_col]]))))

} else if (length(unique(xenium_obj@meta.data[[cluster_col]])) <= 36) {
    cols <- as.vector(pals::polychrome(length(unique(xenium_obj@meta.data[[cluster_col]]))))

} else {
    cols <- as.vector(
        Polychrome::createPalette(
            N = length(unique(xenium_obj@meta.data[[cluster_col]])),
            seedcolors = c("#FF0000", "#00FF00", "#0000FF")
        )
    )
}

# Check if the input was a list of objects or a single object
vln_plot <- VlnPlot(
    xenium_obj,
    layer = "data",
    features = rev(as.character(unique(marker_gene_list))),
    group.by = cluster_col,
    assay = opt$assay,
    stack = TRUE,
    combine = TRUE,
    fill.by = "ident",
    cols = cols,
    flip = TRUE,
    y.max = 1.5,
    same.y.lims = F,
    pt.size = 0.00005
) + NoLegend()

###################
### OUTPUT PLOT ###
###################

# Output the plot
ggsave(
    opt$outfile,
    plot = vln_plot,
    width = opt$width,
    height = opt$height,
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
