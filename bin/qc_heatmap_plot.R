#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package
library(stringr)
library(tidyr)
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
        help="Gene Pair stats to be analyzed"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="heatmap_plot.png",
        metavar="path",
        help="The output name for the heatmap plot"),
    make_option(
        c("--width"),
        type="integer",
        default=2000,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=2000,
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
        c("--ncols_vln"),
        type="integer",
        default=1,
        help="Number of cols for the plot (if there are multiple features)"),
    make_option(
        c("--pt_size"),
        type="double",
        default=0.1,
        help="Size of the points"),
    make_option(
        c("--alpha"),
        type="double",
        default=0.5,
        help="Alpha value for the points")
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
gene_pair_stat_info <- read.csv(
    file = opt$input
)

#####################
### BARNYARD PLOT ###
#####################

# Count occurrences of Gene pairs and FilterPass
pair_counts <- gene_pair_stat_info %>%
    mutate(pair = paste(gene1, gene2, sep = "_")) %>%
    group_by(pair, filter_pass) %>%
    summarise(count = n(), .groups = "drop") %>%
    separate(pair, into = c("gene1", "gene2"), sep = "_")

# Create a matrix for the heatmap
heatmap_data <- pair_counts %>%
    pivot_wider(names_from = filter_pass, values_from = count, values_fill = 0) %>%
    rename(count_true = `TRUE`, count_false = `FALSE`) %>%
    mutate(color = ifelse(count_true > count_false, "red", "blue"))

# Create the heatmap data frame
heatmap_df <- heatmap_data %>%
    pivot_longer(cols = starts_with("count_"), names_to = "status", values_to = "count") %>%
    mutate(status = ifelse(status == "count_true", "True", "False"))

# Plot heatmap using ggplot2
ggplot(heatmap_df, aes(x = gene1, y = gene2, fill = color)) +
    geom_tile() +
    scale_fill_identity() +
    theme_minimal() +
    labs(title = "Gene Pair Heatmap", x = "Gene1", y = "Gene2")

ggsave(filename = opt$outfile, width = opt$width, height = opt$height, units = "px", dpi = 300)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
