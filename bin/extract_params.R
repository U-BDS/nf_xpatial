#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(optparse)          # For parsing commandline arguments
library(Seurat)            # For handling Seurat objects
library(SpatialExperiment) # For handling SpatialExperiment objects
library(stringr)           # For string manipulation

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
        c("-o", "--outfile"),
        type="character",
        default="params.tsv",
        metavar="path",
        help="The extracted params")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
spe_xenium_obj <- readRDS(file = opt$input)

######################
### EXTRACT PARAMS ###
######################

extract_named_list <- function(params_list, source_name = NULL) {
    if (!is.list(params_list) || is.null(names(params_list))) {
        warning("Metadata element is not a named list.")
        return(NULL)
    }

    # Collapse values for each parameter and store in a data frame
    data.frame(
        Parameter = names(params_list),
        Value = sapply(params_list, function(x) paste(x, collapse = ", ")),
        Source = if (!is.null(source_name)) source_name else NA_character_,
        stringsAsFactors = FALSE
    )
}

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

all_params <- lapply(names(spe_xenium_obj@metadata), function(name) {
    extract_named_list(spe_xenium_obj@metadata[[name]], source_name = name)
})

result_df <- unique(do.call(rbind, all_params[!sapply(all_params, is.null)]))

###################
### WRITE TABLE ###
###################

# Save the params
write.table(
    all_params,
    file = opt$outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
