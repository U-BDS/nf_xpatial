#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(arrow)      # For reading/writing Parquet files
library(dplyr)      # For data manipulation
library(jsonlite)   # For working with JSON data
library(knitr)      # For generating reports
library(optparse)   # For parsing commandline arguments
library(progressr)  # For progress bars
library(purrr)      # Functional programming tools
library(Seurat)     # Main analysis package

# Plotting
library(patchwork)  # For combining plots

# Set options
options(future.globals.maxSize = 8192 * 1024^2)

############################
### FUNCTION DEFINITIONS ###
############################

#' process_json_metadata_and_df
#' This function extracts relevant metadata and gene data from a JSON file and flattens nested structures.
#' @param json_file Path to a JSON file.
#'
#' @return metadata Species, tissue, and gene panel details
#' @return df Gene-specific data including gene coverage, category, gene name, and ID.

process_json_metadata_and_df <- function(json_file) {
  # Read the JSON file
  data <- fromJSON(json_file)
  
  # Extract metadata
  metadata <- list(
    Species = data$payload$panel$species,
    Tissue = data$payload$panel$tissue,
    Base_panel_num_gene_targets = as.integer(data$payload$panel$type$data$base_panel_num_gene_targets),
    Num_gene_targets = as.integer(data$payload$panel$num_gene_targets),
    Total_num_gene_targets = as.integer(data$payload$panel$type$data$base_panel_num_gene_targets) +
      as.integer(data$payload$panel$num_gene_targets)
  )
  
  # Normalize JSON data and handle the codewords column
  df <- as.data.frame(data$payload$targets) %>%
    mutate(codewords = sapply(codewords, function(x) if (length(x) == 1) x else NA))
  
  df$gene_coverage <- df$info$gene_coverage
  df$category <- df$source$category
  df$GeneName <- df$type$data$name
  df$id <- df$type$data$id
  df <- df[c("gene_coverage", "category", "GeneName", "id")]
  
  return(list(metadata = metadata, df = df))
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
    warning("'flip.xy' column is missing from metatdata. 'flip.xy' is set to default \n")
    sample_metadata$flip.xy <- NA
}

# Print out columns included in metadata that are not required
addl_metadata_cols <- setdiff(colnames(sample_metadata), req_cols)
print(addl_metadata_cols)

# only grab the row containing the sample
sample_metadata <- sample_metadata[sample_metadata$SampleID == opt$sample, ]

nrow(sample_metadata)
sample_metadata

if (nrow(sample_metadata) > 1) {
    stop(paste("The metadata should only contain one row for each sample. Sample ", opt$sample, "occurs multiple times"))
}

if (nrow(sample_metadata) < 1) {
    stop(paste("The metadata should only contain one row for each sample. Sample ", opt$sample, "is not present"))
}

sample_metadata_row <- sample_metadata[1,]

############################
### LOAD XENIUM METADATA ###
############################

print(paste("Loading ", opt$input, "/cell_feature_matrix"))
# Load xenium object
xenium.obj <- if (!is.na(sample_metadata_row[["flip.xy"]]) && !is.null(sample_metadata_row[["flip.xy"]])) {

    LoadXenium(
        opt$input,
        fov = "fov",
        segmentations = "cell",
        flip.xy = sample_metadata_row[["flip.xy"]]
    )

} else {

    LoadXenium(
        opt$input,
        fov = "fov",
        segmentations = "cell"
    )

}

# Replace underscores with dashes in feature names
rownames(xenium.obj) <- gsub("_", "-", rownames(xenium.obj))

# Remove cells with 0 counts
message("nrow before filtering 0 count cells: ", nrow(xenium.obj@meta.data))
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
message("nrow after filtering 0 count cells: ", nrow(xenium.obj@meta.data))

# Assign project metadata
project <- sample_metadata_row[["SampleID"]]

xenium.obj@project.name <- project
xenium.obj$Sample <- project

xenium.obj$orig.ident <- project
Idents(xenium.obj) <- xenium.obj$orig.ident

# Add additional metadata
for (col in setdiff(names(sample_metadata_row), append(req_cols, "flip.xy"))) {
    xenium.obj[[col]] <- as.character(sample_metadata_row[[col]])
    xenium.obj@meta.data[[col]] <- as.factor(xenium.obj@meta.data[[col]])
}

###############################
### PARSE PANEL AND FEATURE ###
###############################

# Extract gene names
gene_names <- rownames(xenium.obj)

# Process the gene panel json
gene_panel <- process_json_metadata_and_df(
    file.path(paste(opt$input, "/gene_panel.json", sep=""))
)

if (is.null(xenium.obj@misc)) {
    xenium.obj@misc <- list()
}
xenium.obj@misc$panel_metadata <- gene_panel$metadata

print(xenium.obj@misc)

# Process teh features csv
print("Reading in features file")
df_features <- read.csv(
    file.path(paste(opt$input, "/cell_feature_matrix/features.tsv.gz", sep="")),
    sep = "\t",
    header = FALSE,
    col.names = c("Ensembl", "GeneName", "Class")
)

# Filter and reorder features and gene_panel based on the gene_names
print("Filtering gene_panel")
df_gene <- gene_panel$df %>%
    filter(GeneName %in% gene_names) %>% # Filter by gene_names
    arrange(match(GeneName, gene_names)) # Reorder to match gene_name order

print("Filtering features")
df_features <- df_features %>%
    filter(GeneName %in% gene_names) %>% # Filter by gene_names
    arrange(match(GeneName, gene_names)) # Reorder to match gene_name order

# Add combined DataFrame to Xenium object's metadata
xenium.obj@assays$Xenium@meta.data <- left_join(df_gene, df_features, by = "GeneName")

# Display Results
print("Metadata added to Xenium object:")
print(head(xenium.obj@assays$Xenium@meta.data))
print(tail(xenium.obj@assays$Xenium@meta.data))


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

####################
### ADD CELL IDS ###
####################

xenium.obj@meta.data$Cell_ID <- rownames(xenium.obj@meta.data)

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
