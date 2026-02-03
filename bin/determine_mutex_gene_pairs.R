#!/usr/bin/env Rscript

set.seed(1234)

######################
### LOAD LIBRARIES ###
######################

# General Utilities
library(dplyr)      # For data manipulation
library(optparse)   # For parsing commandline arguments
library(purrr)      # Functional programming tools
library(Seurat)     # Main analysis package

###############################
### COMMAND-LINE PARAMETERS ###
###############################

params_list <- list(
    make_option(
        c("-g", "--gene_pair_stats"),
        type="character",
        default=NULL,
        metavar="path",
        help="The gene pair stats file"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="mutex_gene_pairs.csv",
        metavar="path",
        help="The output name for the mutually exclusive gene pairs file"),
    make_option(
        c("--max_spearman_p"),
        type="double",
        default=0.05,
        help="The maximum spearman p value"),
    make_option(
        c("--max_spearman_r"),
        type="double",
        default=-0.2,
        help="The maximum spearman r value"),
    make_option(
        c("--min_gini"),
        type="double",
        default=0.7,
        help="The minimum gini index")
    )

opt_parser <- OptionParser(option_list=params_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene_pair_stats)) {
    print_help(opt_parser)
    stop("Please provide the gene pair stats file as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

gene_pair_stats_df <- read.csv(opt$gene_pair_stats, sep = ",")

############################
### GET MUTEX GENE PAIRS ###
############################

mutex_gene_pair_df <- gene_pair_stats_df %>%
    # Filter by spearman correlation
    filter(
        spearman_p < opt$max_spearman_p,
        spearman_r < opt$max_spearman_r
    ) %>%
    # Filter by Gini index
    filter(
        gini_gene1 > opt$min_gini,
        gini_gene2 > opt$min_gini
    ) %>%
    # Sort by spearman correlation (most negative first)
    arrange(spearman_r)

#################
### SAVE DATA ###
#################

write.csv(
    mutex_gene_pair_df,
    opt$outfile,
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
