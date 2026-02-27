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
        default="vln_plot.pdf",
        metavar="path",
        help="The output name for the image"),
    make_option(
        c("--max_gene_per_group"),
        type="integer",
        default=50,
        help="The maximum number of genes to plot per group")
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

# Read in teh marker list
marker_list <- read.csv(opt$marker_list, sep = ",")

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
    opt$outfile, width = 15, height = 15
)

for (group_name in unique(marker_list$group)) {
    genes <- marker_list[marker_list$group == group_name,]$gene

    # We need to split the gene list to avoid overcrowding the plot
    gene_groups <- split(genes, ceiling(seq_along(genes) / opt$max_gene_per_group))

    vln_plot_title <- paste0(opt$cluster_col, "\n", group_name)

    group_num = 1
    for (gene_group in gene_groups) {

        # Update the title to include the group number if there are multiple groups
        vln_plot_title <- paste0(vln_plot_title, ifelse(length(gene_groups) > 1, paste0(" - ", group_num), ""))

        vln_plot <- VlnPlot(
            xenium_obj,
            layer = "data",
            features = gene_group,
            group.by = opt$cluster_col,
            assay = opt$assay,
            stack = TRUE,
            combine = TRUE,
            fill.by = "ident",
            cols = cols,
            flip = TRUE,
            y.max = 3,
            same.y.lims = F,
            pt.size = 0.00005
        ) + RotatedAxis() + ggtitle(vln_plot_title)

        group_num = group_num + 1
    }

    # Output the plot
    print(vln_plot)
}
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
