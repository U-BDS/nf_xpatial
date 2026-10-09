#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)   # For parsing commandline arguments
library(purrr)      # Functional programming tools
library(Seurat)     # Main analysis package

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium results to be analyzed"),
    make_option(
        c("-m", "--metadata"),
        type="character",
        default=NULL,
        metavar="path",
        help="The metadata file to be incoporated into the seurat object"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="seurat_obj.rds",
        metavar="path",
        help="The output name for the seurat object"),
    make_option(
        c("-s", "--sample"),
        type="character",
        default=NULL,
        metavar="string",
        help="The name of the sample")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}
opt$input <- trimws(opt$input)

if (is.null(opt$metadata)) {
    print_help(opt_parser)
    stop("Please provide the metadata file to add to the seurat object.", call. = FALSE)
}

if (is.null(opt$sample)) {
    print_help(opt_parser)
    stop("Please provide the name of the sample.", call. = FALSE)
}

###########################
### LOAD METADATA INPUT ###
###########################

sample_metadata <- read.csv(opt$metadata)

# TODO: Move to nextflow
# TODO: Figure out what columns are needed
# Validate metadata
# Define required columns
req_cols <- c("SampleID")

# Check for missing columns
missing_cols <- setdiff(req_cols, colnames(sample_metadata))
if (length(missing_cols) > 0) {
    stop(paste("Missing required columns - ", paste(missing_cols, collapse = ", ")))
}

# Check if "flip.xy" is missing, if it is issue a warning and set it to NA
if (!"flip.xy" %in% colnames(sample_metadata)) {
    warning("'flip.xy' column is missing from metadata. 'flip.xy' is set to default \n")
    sample_metadata$flip.xy <- NA
}

# Print out columns included in metadata that are not required
addl_metadata_cols <- setdiff(colnames(sample_metadata), req_cols)
print(addl_metadata_cols)

# only grab the row containing the sample
sample_metadata <- sample_metadata[sample_metadata$SampleID == opt$sample, ]

if (nrow(sample_metadata) > 1) {
    stop(paste("The metadata should only contain one row for each sample. Sample ", opt$sample, "occurs multiple times"))
}

if (nrow(sample_metadata) < 1) {
    stop(paste("The metadata should only contain one row for each sample. Sample ", opt$sample, "is not present"))
}

sample_metadata_row <- sample_metadata[1,]

########################
### LOAD XENIUM DATA ###
########################

# Load xenium object
xenium.obj <- readRDS(file = opt$input)

# Set size limit to object size plus a little
options(future.globals.maxSize = as.numeric(object.size(xenium.obj)) * 3.2)

# Add additional metadata
for (col in setdiff(names(sample_metadata_row), append(req_cols, "flip.xy"))) {
    xenium.obj[[col]] <- as.character(sample_metadata_row[[col]])
    xenium.obj@meta.data[[col]] <- as.factor(xenium.obj@meta.data[[col]])
}
print("Added metdata columns to xenium")
###########################
### ADD SAMPLE METADATA ###
############################

# Replace underscores with dashes in feature names
rownames(xenium.obj) <- gsub("_", "-", rownames(xenium.obj))
print("Replaced underscores with dashes in feature names")

# Assign sample name
project <- opt$sample

xenium.obj@project.name <- project
xenium.obj$Sample <- project

xenium.obj$orig.ident <- project
Idents(xenium.obj) <- xenium.obj$orig.ident
print("Added sample name to various places in metadata")

### ADD CELL IDS
xenium.obj@meta.data$Cell_ID <- rownames(xenium.obj@meta.data)
print("Added cell IDs to xenium object")

### ADD CELL COUNT
Misc(xenium.obj, slot = "cell_count") <- sum(table(xenium.obj@meta.data$orig.ident))
print("Added cell count to xenium object")

#####################
### ADD CELL AREA ###
#####################

polygons <- xenium.obj@images$fov$segmentations@polygons
area_list <- map_dbl(names(polygons), ~ polygons[[.]]@area)
names(area_list) <- names(polygons)

xenium.obj <- AddMetaData(
    xenium.obj,
    metadata = area_list,
    col.name = "Cell_Area"
)

#################
### SAVE DATA ###
#################

saveRDS(
    object = xenium.obj,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
