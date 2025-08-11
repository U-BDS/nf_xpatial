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
        default="contour_dim_plot.png",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("-e", "--embedding"),
        type="character",
        default="umap",
        help="The embedding to use for the plot"),
    make_option(
        c("-p", "--pt_size"),
        type="numeric",
        default=0.1,
        help="The point size for the plot"),
    make_option(
        c("-c", "--metadat_col"),
        type="character",
        default="seurat_clusters",
        help="The metadata column to use"),
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

##########################
### FUNCTION DEFITIONS ###
##########################
Dimplot_contour_ggplot <- function(seurat_object, Embedding="umap", Rug=F, Metadat_column = NULL, Point_Size = 1){
  
  # Extract embeddings (e.g., tSNE or UMAP) from the Seurat object
  embedding <- Embeddings(seurat_object, reduction = Embedding)  # Replace "umap" with "tsne" or the appropriate reduction
  
  # Extract cluster identities
  if (is.null(Metadat_column)) {
    clusters <- Idents(seurat_object)
  } else {
    clusters <- seurat_object@meta.data[[Metadat_column]]
  }
  
  
  # Create a data frame for plotting
  plot_data <- data.frame(embedding, cluster = clusters)
  # print(head(plot_data))
  
  # Create the plot with ggplot2
  if (Embedding=="umap"){
    message("UMAP")
    # Calculate the center of each cluster for labeling
    label_data <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = plot_data, FUN = mean)
    Plot <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = cluster)) +  # Replace UMAP_1 and UMAP_2 with appropriate column names
      geom_point() + geom_density_2d(color = "black", size = 0.5, alpha=.5) + 
      geom_text(data = aggregate(cbind(umap_1, umap_2) ~ cluster, data = plot_data, mean),
                aes(label = cluster), 
                size = Point_Size, 
                color = "black", 
                hjust = 0.5, 
                vjust = 0.5,
                box.padding = 0.3, 
                point.padding = 0.3) + 
      geom_label(data = label_data, aes(label = cluster), color = "black", fill = "white", size = 5, label.padding = unit(0.2, "lines")) # Add labels in boxes
    # stat_ellipse(aes(group = cluster), type = "norm", linetype = 1, size = 0.5, color = "red") +  # Add irregular ellipses
  } 
  if(Embedding=="tsne"){
    message("t-SNE")
    # Calculate the center of each cluster for labeling
    label_data <- aggregate(cbind(tSNE_1, tSNE_2) ~ cluster, data = plot_data, FUN = mean)
    Plot <- ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +  # Replace UMAP_1 and UMAP_2 with appropriate column names
      geom_point() + geom_density_2d(color = "black", size = 0.5, alpha=.5) +
      geom_text(data = aggregate(cbind(tSNE_1, tSNE_2) ~ cluster, data = plot_data, mean),
                aes(label = cluster), 
                size = 5, 
                color = "black", 
                hjust = 0.5, 
                vjust = 0.5,
                box.padding = 0.3, 
                point.padding = 0.3) +
      geom_label(data = label_data, aes(label = cluster), color = "black", fill = "white", size = 5, label.padding = unit(0.2, "lines")) # Add labels in boxes
    
    # stat_ellipse(aes(group = cluster), type = "norm", linetype = 1, size = 0.5, color = "red") +  # Add irregular ellipses
  }
  
  
  if(Rug==TRUE){
    return(Plot + geom_rug() + theme_bw())
  }
  
  return(Plot+theme_bw())
  
}

##################
### LOAD INPUT ###
##################

xenium_obj <- readRDS(
    file = opt$input
)

#################
#### DIM PLOT ###
#################

# TODO: How to do colors
# Need to get all the groups and assign them a color

# Output the plot
dim_plot_contour <- Dimplot_contour_ggplot(
    xenium_obj,
    Embedding = opt$embedding,
    Metadat_column = opt$metadat_col,
    Point_Size = opt$pt_size
)

png(
    opt$outfile,
    width = opt$width,
    height = opt$height
)

plot(dim_plot_contour)
dev.off()

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
