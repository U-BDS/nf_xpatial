#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # Commandline arguments
library(Seurat)     # Main analysis package
library(dplyr)      # Data manipulation

# Plotting
library(ggplot2)
library(patchwork)
library(tools)
library(scales)
library(ggplotify)
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
        default="proportion_plot.png",
        metavar="path",
        help="The output name for the proportion plot"),
    make_option(
        c("-f", "--fill_col"),
        type="character",
        default=NULL,
        help="Name of the column to plot"),
    make_option(
        c("-x", "--x_col"),
        type="character",
        default=NULL,
        help="Name of the column to plot on the x-axis"),
    make_option(
        c("--width"),
        type="integer",
        default=10,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=7,
        help="Height of the plot"),
    make_option(
        c("--dpi"),
        type="integer",
        default=300,
        help="The dpi for the plot")
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

xenium_objs <- readRDS(
    file = opt$input
)

########################
#### PROPORTION PLOT ###
########################

# Combine relevant metadata columns from each object
combined_metadata <- do.call(rbind, lapply(1:length(xenium_objs), function(i) {
    # Get the metadata for the current object
    metadata <- xenium_objs[[i]]@meta.data

    # Add sample (project_name) as column
    metadata$Sample <- xenium_objs[[i]]@project.name

    # select relevant columns
    classification_df <- metadata %>%
        dplyr::select(Cell_ID, !!opt$fill_col, !!opt$x_col)
    
    return(classification_df)
}))

# Plot the data
plot <- ggplot(combined_metadata, aes_string(x = opt$x_col, fill = opt$fill_col)) +
  geom_bar(position = "fill", stat = "count") + # Stacked bar plot with proportions
  labs(title = paste(tools::toTitleCase(gsub("_", " ", opt$fill_col)), "Distribution Across Samples"),
       x = "Sample",
       y = "Proportion",
       fill = tools::toTitleCase(gsub("_", " ", opt$fill_col))) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.direction = "vertical") +
  scale_y_continuous(labels = scales::percent)

# Save the plot
ggsave(
    filename = opt$outfile,
    plot = plot,
    width = opt$width,
    height = opt$width,
    dpi = opt$dpi
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
