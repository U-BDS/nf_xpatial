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
        c("-c", "--cluster_col"),
        type="character",
        default="seurat_clusters",
        help="The column name to pull for cluster numbers"),
    make_option(
        c("-p", "--pt_size"),
        type="numeric",
        default=0.1,
        help="The point size for the plot"),
    make_option(
        c("-s", "--cluster_label_size"),
        type="numeric",
        default=12,
        help="The size of the cluster label for the plot"),
    make_option(
        c("-a", "--axis_label_size"),
        type="numeric",
        default=30,
        help="The size of the axis label for the plot"),
    make_option(
        c("-l", "--legend_label_size"),
        type="numeric",
        default=30,
        help="The size of the legend label for the plot"),
    make_option(
        c("-g", "--element_guide_size"),
        type="numeric",
        default=15,
        help="The size of the legend element for the plot"),    
    make_option(
        c("--width"),
        type="integer",
        default=10000,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=10000,
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
Dimplot_contour_ggplot <- function(seurat_object, Embedding="umap", Rug=F, Cluster_column = NULL, Point_Size = 0.01){
  
  # Extract embeddings (e.g., tSNE or UMAP) from the Seurat object
  embedding <- Embeddings(seurat_object, reduction = Embedding)  # Replace "umap" with "tsne" or the appropriate reduction
  
  # Extract cluster identities
  if (is.null(Cluster_column)) {
    clusters <- Idents(seurat_object)
  } else {
    clusters <- seurat_object@meta.data[[Cluster_column]]
  }
  
  
  # Create a data frame for plotting
  plot_data <- data.frame(embedding, cluster = clusters)
  # print(head(plot_data))
  embedding_1 <- paste0(gsub("_", "", Embedding), "_1")
  embedding_2 <- paste0(gsub("_", "", Embedding), "_2")
  
  # Create the plot with ggplot2
  if (grepl("umap", tolower(Embedding))){
    message("UMAP")
    # Calculate the center of each cluster for labeling
    #label_data <- aggregate(cbind(umap_1, umap_2) ~ cluster, data = plot_data, FUN = mean)
    label_data <- aggregate(
      plot_data[, c(embedding_1, embedding_2)],
      by = list(cluster = plot_data$cluster),
      FUN = mean
    )

    Plot <- ggplot(plot_data, aes(x = .data[[embedding_1]], y = .data[[embedding_2]], color = cluster)) +  # Replace UMAP_1 and UMAP_2 with appropriate column names
      geom_point() + geom_density_2d(color = "black", linewidth = 0.5, alpha=.5) + 
      geom_text(data = label_data,
                aes(label = cluster), 
                size = Point_Size, 
                color = "black", 
                hjust = 0.5, 
                vjust = 0.5) + 
      geom_label(data = label_data, aes(label = cluster), color = "black", fill = "white", size = opt$cluster_label_size, label.padding = unit(0.2, "lines")) # Add labels in boxes
      # stat_ellipse(aes(group = cluster), type = "norm", linetype = 1, size = 0.5, color = "red") +  # Add irregular ellipses
  } 
  if(Embedding=="tsne"){
    message("t-SNE")
    # Calculate the center of each cluster for labeling
    label_data <- aggregate(cbind(tSNE_1, tSNE_2) ~ cluster, data = plot_data, FUN = mean)
    Plot <- ggplot(plot_data, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +  # Replace UMAP_1 and UMAP_2 with appropriate column names
      geom_point() + geom_density_2d(color = "black", linewidth = 0.5, alpha=.5) +
      geom_text(data = aggregate(cbind(tSNE_1, tSNE_2) ~ cluster, data = plot_data, mean),
                aes(label = cluster), 
                size = 5, 
                color = "black", 
                hjust = 0.5, 
                vjust = 0.5) +
      geom_label(data = label_data, aes(label = cluster), color = "black", fill = "white", size = opt$cluster_label_size, label.padding = unit(0.2, "lines")) # Add labels in boxes
    
    # stat_ellipse(aes(group = cluster), type = "norm", linetype = 1, size = 0.5, color = "red") +  # Add irregular ellipses
  }
  
  
  if(Rug==TRUE){
    return(Plot + 
             geom_rug() + 
             theme_bw() +
             theme(legend.title = element_text(size = opt$legend_label_size),
                   legend.text = element_text(size = opt$legend_label_size),
                   axis.title = element_text(size = opt$axis_label_size),
                   axis.text = element_text(size = opt$axis_label_size)) +
             guides(color = guide_legend(override.aes = list(size = opt$element_guide_size)))
           ) 
  }
  
  return(Plot +
           theme_bw() +
           theme(legend.title = element_text(size = opt$legend_label_size),
                 legend.text = element_text(size = opt$legend_label_size),
                 axis.title = element_text(size = opt$axis_label_size),
                 axis.text = element_text(size = opt$axis_label_size)) +
           guides(color = guide_legend(override.aes = list(size = opt$element_guide_size)))
         )
  
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
    Cluster_column = opt$cluster_col,
    Point_Size = opt$pt_size
) + plot_annotation(
        opt$cluster_col,
        theme = theme(
            plot.title = element_text(size = 32, hjust = 0.5, face = "bold")
        )) + theme(plot.title = element_text(size = opt$axis_label_size, hjust = 0.5))

###################
### OUTPUT PLOT ###
###################

# Output the plot
ggsave(
    opt$outfile,
    plot = dim_plot_contour,
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
