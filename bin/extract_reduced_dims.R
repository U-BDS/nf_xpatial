#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)          # For parsing commandline arguments
library(Seurat)            # For handling Seurat objects
library(stringr)           # For string manipulation
library(SpatialExperiment) # For handling SpatialExperiment objects

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium object"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to use"),
    make_option(
        c("-p", "--param_string"),
        type="character",
        help="The parameter string to use to identify clustering run"),
    make_option(
        c("-c", "--clustering_method"),
        type="character",
        help="The clustering method used (e.g., BANKSY, Harmony)"),
    make_option(
        c("-f", "--filter"),
        type="character",
        help="Comma delimited list of dimensions to skip"
    ),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="seurat_reduced_dims.csv",
        metavar="path",
        help="The suffix to be used for csv's containing seurat reduced dimensions")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

############################
### EXTRACT REDUCED DIMS ###
############################

if (class(xenium_obj) == "Seurat") {
    DefaultAssay(xenium_obj) <- opt$assay

    if (length(Reductions(xenium_obj)) == 0) {
        stop("No reductions found in xenium object.")
    }

    reductions <- Reductions(xenium_obj)
    # Filter out reductions if specified
    if (!is.null(opt$filter)) {
        filter_dims <- str_split_1(opt$filter, ",")

        reductions <- reductions[!reductions %in% filter_dims]
    }

    # Save each reduction-specific data frame to a csv file
    for (reduc in reductions) {
        
        file_prefix <- paste0(gsub("_","-",reduc), "_", opt$param_string)

        # Extract embeddings
        extracted_embeddings <- Embeddings(xenium_obj, reduction = reduc)

        embeddings_df <- data.frame(
            Index = rownames(extracted_embeddings),
            extracted_embeddings,
            stringsAsFactors = FALSE
        )

        write.csv(
            embeddings_df,
            file = paste0(file_prefix, "_embeddings_", opt$outfile),
            row.names = FALSE,
            quote = FALSE
        )

        # NOTE: Banksy objects will not have embeddings or loadings

        if (opt$clustering_method != "BANKSY") {
            # Extract Loadings
            extracted_loadings <- Loadings(xenium_obj, reduction = reduc)

            loadings_df <- data.frame(
                Index = rownames(extracted_loadings),
                extracted_loadings,
                stringsAsFactors = FALSE
            )

            write.csv(
                loadings_df,
                file = paste0(file_prefix, "_loadings_", opt$outfile),
                row.names = FALSE,
                quote = FALSE
            )

            # Extract Stdev
            extracted_stdev <- Stdev(xenium_obj, reduction = reduc)

            stdev_df <- data.frame(
                extracted_stdev,
                stringsAsFactors = FALSE
            )

            write.csv(
                stdev_df,
                file = paste0(file_prefix, "_stdev_", opt$outfile),
                row.names = FALSE,
                quote = FALSE
            )
        }
    }

} else if (class(xenium_obj) == "SpatialExperiment") {
    # Set default assay
    mainExpName(xenium_obj) <- opt$assay

    reductions <- reducedDimNames(xenium_obj)
    # Filter out reductions if specified
    if (!is.null(opt$filter)) {
        filter_dims <- str_split_1(opt$filter, ",")

        reductions <- reductions[!reductions %in% filter_dims]
    }

    # Save each sample-specific data frame to a csv file
    for (reduc in reductions) {
        file_prefix <- paste0(gsub("_","-",reduc), "_", opt$param_string)
        reduced_dim_data <- reducedDim(xenium_obj, reduc)

        if (reduc == "UMAP_BANKSY_harmony") {
            colnames(reduced_dim_data) <- c("umap_1", "umap_2")
        } 

        reduced_dim_df <- data.frame(
            Index = rownames(reduced_dim_data),
            reduced_dim_data,
            stringsAsFactors = FALSE
        )

        write.csv(
            reduced_dim_df,
            file = paste0(file_prefix, "_embeddings_", opt$outfile),
            row.names = FALSE
        )
    }

} else {
    stop("Input xenium object is not of class SeuratObject or SpatialExperiment.")
}

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
