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
        c("-b", "--banksy_clust_info"),
        type="character",
        default=NULL,
        metavar="path",
        help="The cluster information file"),
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

if (is.null(opt$banksy_clust_info)) {
    print_help(opt_parser)
    stop("Please provide the banksy cluster file to add to the xenium object.", call. = FALSE)
}

##################
### LOAD INPUT ###
##################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Read in the banksy cluster info
metadata_df <- as.data.frame(read.table(
    file = opt$banksy_clust_info,
    header = TRUE,
    sep = ","
))

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

    # 1) Make the chosen column a factor & set as Idents
    colnames(seurat_object@meta.data)
    seurat_object@meta.data[[cname]] <- as.factor(seurat_object@meta.data[[cname]])
    Idents(seurat_object) <- cname

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

    # 6) Plot UMAP
    seurat_object_temp <- seurat_object
    RUN_ID <- cname
    colnames(clust_df)[2:3] <- c("UMAP_1", "UMAP_2")

    # clust_df
    Rownames <- clust_df$Index

    # # Ensure that 'Index' column is set as row names in the UMAP dataframe
    rownames(clust_df) <- clust_df$Index
    clust_df <- as.matrix(clust_df[, -1])  # Remove the 'Index' column to keep only UMAP coordinates
    rownames(clust_df) <- Rownames
    # 
    # Now we add the external UMAP data as a new reduction in the Seurat object
    clust_df
    assay_name
    seurat_object_temp[[RUN_ID]] <- CreateDimReducObject(
        embeddings = clust_df,  # UMAP coordinates as a matrix
        key = "UMAP_",  assay = assay_name# Prefix for the dimension names (optional, e.g., UMAP_1, UMAP_2)
    )

    # To verify if the new reduction has been added correctly
    # seurat_object_temp[[RUN_ID]]
    UMAP_Plot <- DimPlot(
        seurat_object_temp,
        reduction = RUN_ID,
        pt.size = 0.1,
        split.by = "Sample",
        cols = colors ,
        combine = F,
        raster = F,
        ncol=2
    )
    UMAP_Plot <- UMAP_Plot[[1]] + ggplot2::theme(legend.position = "none")
    length(UMAP_Plot)

    rm(seurat_object_temp)
    gc()

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

#######################
#### QC_BANKSY PLOT ###
#######################

cnames <- colnames(xenium_obj@meta.data)
cnames <- cnames[grep("^clust", cnames)]

for (cname in cnames) {
    Plot <- generate_seurat_plots(
        cname = cname,
        seurat_object = xenium_obj,  # Your Seurat object
        clust_df = metadata_df,  # Named list of UMAP CSV files
        output_file = paste0(cname,".",opt$outfile),  # Directory to save plots
        return_plot = FALSE,  # Set to TRUE to return the combined plot
        assay_name = opt$assay,
        gene_list = opt$gene_list
    )
}

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
