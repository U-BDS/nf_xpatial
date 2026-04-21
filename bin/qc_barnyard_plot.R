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

###########################
### FUNCTION DEFINTIONS ###
###########################

# From ubdsR
plot_barnyard <- function(
    xenium_object, genes_of_interest, x_lab, y_lab, 
    colors = c("#000000", "#0072B2", "#D55E00", "#CC79A7"), 
    ignore_none = TRUE, legend = TRUE, 
    plot_types = "density", layer = "data", 
    annotate_plot = FALSE) {
  
    # Plot types options
    plot_types_options = c("density", "histogram", "boxplot", "violin", "densigram")
    
    # Validate Plot_types
    if (!all(plot_types %in% plot_types_options)) {
        plot_types <- NULL
        message("Setting Plot_types to NULL, will plot only scatter plot!")
        message(c("density", "histogram", "boxplot", "violin", "densigram"))
    }

    # Fetch expression data for the genes of interest
    expression_data <- FetchData(xenium_object, vars = genes_of_interest, layer = layer)
    
    # Create a data frame for plotting
    plot_data <- data.frame(
        x_gene = expression_data[[genes_of_interest[1]]],
        y_gene = expression_data[[genes_of_interest[2]]]
    )

    # Classify cells based on gene expression
    plot_data$Category <- "None"
    plot_data$Category[plot_data$x_gene > 0 & plot_data$y_gene == 0] <- x_lab
    plot_data$Category[plot_data$x_gene == 0 & plot_data$y_gene > 0] <- y_lab
    plot_data$Category[plot_data$x_gene > 0 & plot_data$y_gene > 0] <- "Mixed"
    
    if (ignore_none) {
        plot_data <- plot_data[plot_data$Category != "None", ]
    }
    
    # Color mapping
    color_mapping <- setNames(colors, c("None", x_lab, y_lab, "Mixed"))
    
    # Calculate statistics
    contingency_table <- table(plot_data$Category)
    
    # Print the number of cells in each category
    print("Number of cells in each category:")
    print(as.data.frame(contingency_table))  # Convert to data frame for readability
    
    # Convert contingency table to data frame for annotation
    total_cells <- sum(contingency_table)
    annotation_data <- as.data.frame(contingency_table)
    colnames(annotation_data) <- c("Category", "Count")
    annotation_data$Percentage <- round((annotation_data$Count / total_cells) * 100, 2)

    # Create subtitle text
    subtitle_text <- paste(
        paste0(annotation_data$Category, ": ", 
            annotation_data$Count, " (", 
            annotation_data$Percentage, "%)"),
        collapse = "\n"
    )

    # Create the scatter plot
    scatter_plot <- ggplot(plot_data, aes(x = x_gene, y = y_gene, color = Category)) +
        geom_point(size = 1, alpha = 0.5) +
        scale_color_manual(values = color_mapping) +
        theme_minimal() +
        labs(
        x = x_lab, 
        y = y_lab,
        subtitle = subtitle_text  # Add subtitle here
        ) +
        theme(legend.position = "bottom") +
        ggtitle(paste0("Barnyard Plot for ", x_lab, " and ", y_lab))
    
    if (!legend) {
        scatter_plot <- scatter_plot + theme(legend.position = "none")
    }
    
    # Add annotations in the plot if annotate_plot is TRUE
    if (annotate_plot) {
        scatter_plot <- scatter_plot +
        annotate(
            "text",
            x = max(plot_data$x_gene) * 0.8, 
            y = max(plot_data$y_gene) * seq(0.9, 0.5, length.out = nrow(annotation_data)),
            label = paste0(
            annotation_data$Category, ": ", 
            annotation_data$Count, " (", annotation_data$Percentage, "%)"
            ),
            hjust = 0, size = 4, color = "black"
        ) +
        ggtitle(paste0("Barnyard Plot for ", x_gene, " and ", y_gene))
    }

    scatter_plot <- scatter_plot + theme_classic()
    
    # Add marginal plot types if specified
    if (is.null(plot_types)) {
        return(ggplotify::as.ggplot(scatter_plot))
    } else {
        plot <- ggExtra::ggMarginal(scatter_plot, groupColour = TRUE, groupFill = TRUE, type = plot_types)
        return(ggplotify::as.ggplot(plot))
    }
}

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
        default="barnyard_plot.png",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("--gene_pair_stats"),
        type="character",
        default=NULL,
        metavar="path",
        help="Gene pair CSV to be analyzed"),
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
xenium_objs <- readRDS(
    file = opt$input
)

# Read in gene pair stats if provided
gene_pair_stats_df <- read.csv(opt$gene_pair_stats, sep = ",")

#####################
### BARNYARD PLOT ###
#####################

# TODO: How to do colors
# Need to get all the groups and assign them a color

apply(gene_pair_stats_df, 1, function(gene_pair_stat) {
    gene1 <- gene_pair_stat[['gene1']]
    gene2 <- gene_pair_stat[['gene2']]

    barnyard_plot <- plot_barnyard(
        xenium_object = xenium_objs,
        genes_of_interest = c(gene1, gene2),
        x_lab = gene1,
        y_lab = gene2,
        colors = c("#000000", "#0072B2", "#D55E00", "#CC79A7"), # color blind friendly choices
        ignore_none = TRUE,
        legend = FALSE,
        plot_types = "density",
        layer = "counts",
        annotate_plot = FALSE
    )

    ggsave(
        paste0(gene1, "_", gene2, ".", opt$outfile),
        plot = barnyard_plot,
        width = opt$width,
        height = opt$height,
        units = "px",
        bg = "white"
    )
})

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
