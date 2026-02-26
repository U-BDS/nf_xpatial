#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package
library(stringr)
library(dplyr)

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
        c("-m", "--marker_list"),
        type="character",
        default=NULL,
        metavar="path",
        help="Marker list"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to operate on"),
    make_option(
        c("-c", "--cluster_col"),
        type="character",
        default="seurat_clusters",
        help="The column name to pull for cluster numbers"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="split_cluster_plot.png",
        metavar="path",
        help="The output name for the image")
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

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

marker_list <- read.csv(opt$marker_list, sep = "\t")

marker_list <- marker_list[c("Cell_Type", "Subtype", "Marker_Genes")]
#marker_list$Cell_Type <- gsub("/", "", marker_list$Cell_Type)
#marker_list$Cell_Type <- gsub(" ", "_", marker_list$Cell_Type)

#marker_list$Subtype <- gsub("/||\\(||\\)", "", marker_list$Subtype)
#marker_list$Subtype <- gsub("\\+", "_pos", marker_list$Subtype)
#marker_list$Subtype <- gsub(" ", "_", marker_list$Subtype)

marker_list$Marker_Genes <- gsub(" ", "", marker_list$Marker_Genes)

length(marker_list$Subtype)

marker_list[1,]$Cell_Type

#######################
#### QC_BANKSY PLOT ###
#######################

# Get colors
cols <- as.vector(
    Polychrome::createPalette(
    N = length(unique(xenium_obj@meta.data[[opt$cluster_col]])),
    seedcolors = c("#FF0000", "#00FF00", "#0000FF")
    )
)

pdf(
    paste0(opt$cluster_col,".dot_plots.pdf"), width = 15, height = 7
)

for (k in 1:length(marker_list$Cell_Type)) {
    cell_type <- marker_list[k,]$Cell_Type
    subtype <- marker_list[k,]$Subtype
    genes <- str_split_1(marker_list[k,]$Marker_Genes, ",")

    # Check if the input was a list of objects or a single object
    dot_plot <- DotPlot(
        xenium_obj,
        features = genes,
        group.by = opt$cluster_col,
        assay = opt$assay,
        cols = c("red", "blue")
    ) + RotatedAxis() + ggtitle(paste0(opt$cluster_col, "\n", cell_type, " - ", subtype))

    # Output the plot
    print(dot_plot)
}
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
