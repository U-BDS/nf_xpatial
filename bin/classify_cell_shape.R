#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(Seurat)     # Main analysis package
library(geometry)   # For geometric calculations

############################
### FUNCTION DEFINITIONS ###
############################

classify_polygon_shape <- function(xenium_obj, cell_id = NULL) {
    # Get the coordinates of the cell
    coords <- xenium_obj@images$fov$segmentations@polygons[[cell_id]]@Polygons[[1]]@coords
    centroid <- xenium_obj@images$fov$segmentations@polygons[[cell_id]]@Polygons[[1]]@labpt
    
    # Get x and y coordinates from coords
    x_coords <- coords[, 1]
    y_coords <- coords[, 2]

    # Calculate the convex hull of the polygon
    hull <- chull(x_coords, y_coords)
    hull_coords <- data.frame(x = x_coords[hull], y = y_coords[hull])

    # Calculate the area of the polygon
    polygon_area <- abs(sum(x_coords * c(y_coords[-1], y_coords[1]) - y_coords * c(x_coords[-1], x_coords[1])) / 2)
    
    # Calculate the perimeter of the polygon
    perimeter <- sum(sqrt(diff(c(x_coords, x_coords[1]))^2 + diff(c(y_coords, y_coords[1]))^2))

    # Calculate the aspect ratio
    aspect_ratio <- (max(x_coords) - min(x_coords)) / (max(y_coords) - min(y_coords))
    
    # Calculate the circularity
    circularity <- (4 * pi * polygon_area) / (perimeter^2)
    
    classification <- "Polygonal"
    # Classify the shape based on circularity
    if (circularity > 0.8) {
        classification <- "Circular"

    } else if (circularity > 0.8 && aspect_ratio < 1.05) {
        classification <- "Polygonal"

    } else if (aspect_ratio >= 1.05 && aspect_ratio <= 1.2) {
        classification <- "Polygonal"

    } else if (aspect_ratio > 1.2) {
        classification <- "Elongated"

    }

    return(classification)
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
        help="The input xenium object"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="cell_shape_xenium_obj.rds",
        metavar="path",
        help="The output xenium object"),
    make_option(
        c("-s", "--shape_file"),
        type="character",
        default="cell_shape_classification.csv",
        metavar="path",
        help="The output cell classification file")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the xenium ojbect file as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

###########################
### CLASSIFY CELL SHAPE ###
###########################

cell_ids <- colnames(xenium_obj)

# Iterate through each cell and classify its shape
shape_classification <- sapply(X = seq_along(cell_ids), FUN = function(x) {
  
  classify_polygon_shape(xenium_obj, cell_id = cell_ids[x])
  
})

cell_classification <- data.frame(
    cell_id = cell_ids,
    classification = shape_classification,
    stringsAsFactors = FALSE
)

# Add the classification results to the xenium object
xenium_obj <- AddMetaData(
    object = xenium_obj,
    metadata = cell_classification$classification,
    col.name = "shape_classification"
)
xenium_obj@meta.data

#################
### SAVE DATA ###
#################

# Save the xenium object with the classification
saveRDS(
    object = xenium_obj,
    file = opt$outfile 
)

# Save the classification results to a CSV file
write.csv(
    x = cell_classification,
    file = opt$shape_file,
    row.names = FALSE,
    quote = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
