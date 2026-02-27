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
        default="dot_plot.pdf",
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

pdf(
    opt$outfile, width = 15, height = 7
)

for (group_name in unique(marker_list$group)) {
    genes <- marker_list[marker_list$group == group_name,]$gene

    # We need to split the gene list to avoid overcrowding the plot
    gene_groups <- split(genes, ceiling(seq_along(genes) / opt$max_gene_per_group))

    dot_plot_title <- paste0(opt$cluster_col, "\n", group_name)

    group_num = 1
    for (gene_group in gene_groups) {

        # Update the title to include the group number if there are multiple groups
        dot_plot_title <- paste0(dot_plot_title, ifelse(length(gene_groups) > 1, paste0(" - ", group_num), ""))

        # Check if the input was a list of objects or a single object
        dot_plot <- DotPlot(
            xenium_obj,
            features = genes,
            group.by = opt$cluster_col,
            assay = opt$assay,
            cols = c("red", "blue")
        ) + RotatedAxis() + ggtitle(dot_plot_title)

        # Output the plot
        print(dot_plot)

        group_num = group_num + 1
    }
}
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
