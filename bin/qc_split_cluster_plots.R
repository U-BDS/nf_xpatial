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
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to operate on"),
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

#####################
### FUNCTION DEFS ###
#####################

plotCountsProportionsSingle <- function(input_obj,
                                        metadata_variable,
                                        plot_type = c("proportion", "count"),
                                        value_show_type = c("proportion", "count", "none"),
                                        label_digits = 2,
                                        fill_colors = NULL,
                                        show_legend = TRUE,
                                        show_axis_labels = TRUE,
                                        keep_plot_clean = FALSE,
                                        orientation = c("vertical", "horizontal")) {
  
    plot_type <- match.arg(plot_type)
    value_show_type <- match.arg(value_show_type)
    orientation <- match.arg(orientation)

    # Step 1: Fetch metadata and count each group
    metadata <- FetchData(input_obj, vars = metadata_variable) %>%
        mutate(across(everything(), as.factor)) %>%
        group_by(.data[[metadata_variable]]) %>%
        summarize(total = n(), .groups = "drop")

    # Step 2: Compute total and proportions
    total_all <- sum(metadata$total)
    metadata <- metadata %>%
        mutate(proportion = total / total_all)

    # Step 3: Set y-axis value
    if (plot_type == "proportion") {
        y_var <- "proportion"
        y_label <- "Proportion"
    } else {
        y_var <- "total"
        y_label <- "Cell/Spot Count"
    }

    # Step 4: Set label values
    if (value_show_type == "proportion") {
        label_value <- round(metadata$proportion, digits = label_digits)
    } else if (value_show_type == "count") {
        label_value <- metadata$total
    } else {
        label_value <- NULL
    }

    # Step 5: Build base plot (always use vertical mapping)
    plot <- ggplot(metadata, aes(x = "", y = .data[[y_var]], fill = .data[[metadata_variable]])) +
        geom_col(width = 0.5, color = "black")

    # Step 6: Add text labels if needed
    if (!is.null(label_value)) {
        text_layer <- geom_text(
            aes(label = label_value),
            position = position_stack(vjust = 0.5),
            size = 3,
            angle = ifelse(orientation == "horizontal", 90, 0))  # rotate if needed
        plot <- plot + text_layer
    }


    # Step 7: Apply custom fill colors
    if (!is.null(fill_colors)) {
        plot <- plot + scale_fill_manual(values = fill_colors)
    }

    # Step 8: Theme adjustments
    if (keep_plot_clean) {
        plot <- plot +
            labs(x = NULL, y = NULL, fill = NULL) +
            theme_void() +
            theme(legend.position = "none")
    } else {
        # Axis labels
        if (show_axis_labels) {
            plot <- plot +
                labs(x = NULL, y = y_label, fill = metadata_variable)
        } else {
            plot <- plot +
                labs(x = NULL, y = NULL, fill = NULL)
        }
        # Remove clutter
        plot <- plot +
            theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank())

        # Legend control
        if (!show_legend) {
            plot <- plot + theme(legend.position = "none")
        }
    }

    # Step 9: Flip orientation if needed
    if (orientation == "horizontal") {
        plot <- plot + coord_flip()
    }

    return(plot)
}

generate_seurat_plots <- function(cname, seurat_object, clust_df, output_file, assay_name, gene_list, return_plot = FALSE) {

    # 2) Define a color palette, ensuring reproducibility
    set.seed(1234)
    colors <- createPalette(
        length(levels(seurat_object@meta.data[, cname])),
        c("#00ffff", "#ff00ff", "#ffff00"),
        prefix = "",
        M = 1000
    )

    names(colors) <- levels(seurat_object@meta.data[[cname]])

    # Create single horizontal stacked bar showing composition
    composition_bar <- plotCountsProportionsSingle(
        input_obj = seurat_object,
        metadata_variable = cname,
        value_show_type = "count",
        keep_plot_clean = TRUE,
        orientation = "horizontal",
        fill_colors = colors  # match to other plots
    )

    # 4) Create spatial dimension plots (returns a list of plots)
    ImageDim_Plot <- ImageDimPlot(
        seurat_object,
        fov = Images(seurat_object),
        dark.background = FALSE,
        group.by = cname,
        cols = colors,
        combine = F,
        size = 0.3
    )

    # 5) Customize each plot
    fovs <- Images(seurat_object)
    fovs <- gsub("fov.", "", fovs)

    for (i in seq_along(ImageDim_Plot)) {
        ImageDim_Plot[[i]] <- ImageDim_Plot[[i]] +
            ggplot2::ggtitle(fovs[i]) +
            ggplot2::theme(legend.position = "none") +
            coord_flip() +
            scale_x_reverse()
    }

    # To verify if the new reduction has been added correctly
    # seurat_object_temp[[RUN_ID]]
    UMAP_Plot <- DimPlot(
        seurat_object,
        reduction = "umap",
        pt.size = 0.1,
        split.by = "Sample",
        cols = colors ,
        combine = F,
        raster = F,
        ncol=2
    )
    UMAP_Plot <- UMAP_Plot[[1]] + ggplot2::theme(legend.position = "none")
    length(UMAP_Plot)

    # 7) Combine the first 8 dimension plots into a grid
    left_col <- wrap_plots(
        ImageDim_Plot,
        ncol = 2    # 2 columns
    )

    # 8. Wrap the 9th plot (the violin) as a separate patch
    right_col <- wrap_plots(
        UMAP_Plot #ImageDim_Plot[[9]]
    )

    # Combine spatial plots (left), UMAP (right)
    top_panel <- wrap_plots(
        left_col,
        right_col,
        ncol = 2,
        widths = c(5, 5)
    )

    # 9. Now arrange these two "wrapped" patches side by side
    # Final full panel: top plots + composition bar below
    combined_plot <- wrap_plots(
        top_panel,
        composition_bar,
        ncol = 1,
        heights = c(12, 1.2)  # adjust as needed
    )

    # 10. Save the combined plot
    ggsave(
        filename = file.path(output_file),
        plot = combined_plot,
        width = 22,       # adjust as needed
        height = 18       # adjust as needed
    )
    # Optionally return the plot
    if (return_plot) {
        return(combined_plot)
    }
}

##################
### LOAD INPUT ###
##################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

#######################
#### QC_BANKSY PLOT ###
#######################

Plot <- generate_seurat_plots(
    cname = "seurat_clusters",
    seurat_object = xenium_obj,  # Your Seurat object
    output_file = opt$outfile,  # Directory to save plots
    return_plot = FALSE,  # Set to TRUE to return the combined plot
    assay_name = opt$assay,
    gene_list = opt$gene_list
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
