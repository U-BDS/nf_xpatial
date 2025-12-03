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
        help="The name of the sample"),
    make_option(
        c("--flip_xy"),
        action="store_true",
        default=FALSE,
        help="Whether to flip the x and y coordinates for the xenium outputs"),
    make_option(
        c("--mols_qv_threshold"),
        type="double",
        default=0,
        help="The minimum QV for transcript molecules to be retains")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium results as input.", call. = FALSE)
}
opt$input <- trimws(opt$input)

if (is.null(opt$sample)) {
    print_help(opt_parser)
    stop("Please provide the name of the sample.", call. = FALSE)
}

########################
### LOAD XENIUM DATA ###
########################

print(paste("Loading ", opt$input, "/cell_feature_matrix"))
# Load xenium object
xenium.obj <- LoadXenium(
        opt$input,
        fov = "fov",
        segmentations = "cell",
        flip.xy = opt$flip_xy,
        mols.qv.threshold = opt$mols_qv_threshold
)

# Replace underscores with dashes in feature names
rownames(xenium.obj) <- gsub("_", "-", rownames(xenium.obj))

# Assign sample name
project <- opt$sample

xenium.obj@project.name <- project
xenium.obj$Sample <- project

xenium.obj$orig.ident <- project
Idents(xenium.obj) <- xenium.obj$orig.ident

####################
### ADD CELL IDS ###
####################

xenium.obj@meta.data$Cell_ID <- rownames(xenium.obj@meta.data)

######################
### ADD CELL COUNT ###
######################
Misc(xenium.obj, slot = "cell_count") <- sum(table(xenium.obj@meta.data$orig.ident))

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

message("--- Printing misc information ---")
print(xenium.obj@misc)

# Process the features csv
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
