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
library(data.table)        # For data table manipulation
library(ggplot2)           # For plotting

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-i", "--input"),
        type="character",
        default=NULL,
        metavar="path",
        help="The xenium object to be filtered"),
    make_option(
        c("-a", "--assay"),
        type="character",
        default=NULL,
        help="The assay to keep during conversion"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="filtered_xenium_obj.rds",
        metavar="path",
        help="The filtered xenium object"),
    make_option(
        c("-s", "--staggered_plot"),
        type="character",
        default="staggered_plot.png",
        metavar="path",
        help="The output name for the image"),
    make_option(
        c("--width"),
        type="integer",
        default=3500,
        help="Width of the plot"),
    make_option(
        c("--height"),
        type="integer",
        default=2000,
        help="Height of the plot")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
spe_xenium_obj <- readRDS(file = opt$input)

######################
### CONVERT TO SPE ###
######################

# Set the default assay on the spatial experiment object
mainExpName(spe_xenium_obj) <- opt$assay

# Extract spatial coordinates
locs_dt <- data.table(
    cbind(
        spatialCoords(spe_xenium_obj),
        sample_id = factor(spe_xenium_obj$orig.ident)
    )
)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")

locs_dt[, sdimx := sdimx - min(sdimx), by = group]
locs_dt[, sdimx := sdimx + group * (max(locs_dt$sdimx) * 1.5)]

locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe_xenium_obj)
spatialCoords(spe_xenium_obj) <- locs

######################
### STAGGERED PLOT ###
######################

# a good sanity check
metadata_staggered <- merge(colData(spe_xenium_obj),
                            spatialCoords(spe_xenium_obj),
                            by = 0)

staggered_plot <- ggplot(metadata_staggered,
                         aes(orig.ident, sdimx)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(opt$assay)

# Output the plot
ggsave(
  opt$staggered_plot,
  plot = staggered_plot,
  width = opt$width,
  height = opt$height,
  units = "px"
)

#################
### SAVE DATA ###
#################

saveRDS(
    object = spe_xenium_obj,
    file = opt$outfile 
)

####################
### SESSION INFO ###
####################

sessioninfo <- "R_sessionInfo.log"

sink(sessioninfo)
sessionInfo()
sink()
