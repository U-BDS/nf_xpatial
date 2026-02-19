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

############################
### FUNCTION DEFINITIONS ###
############################

calculate_gini_score <- function(x) {
    x <- sort(x)
    n <- length(x)
    G <- sum((2 * (1:n) - n - 1) * x) / (n * sum(x))
    return(G)
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
        help="The xenium rds object"),
    make_option(
        c("-g", "--gene_list"),
        type="character",
        default=NULL,
        metavar="path",
        help="The marker gene list"),
    make_option(
        c("-o", "--outfile"),
        type="character",
        default="gene_pair_stats.csv",
        metavar="path",
        help="The output name for the gene pair stats"),
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

if (is.null(opt$input)) {
    print_help(opt_parser)
    stop("Please provide the Xenium rds object as input.", call. = FALSE)
}

###################
### LOAD INPUTS ###
###################

# Read in xenium_obj
xenium_obj <- readRDS(file = opt$input)

# Read in marker gene list if passed in
marker_gene_list <- NULL
if (!is.null(opt$gene_list) && opt$gene_list != "") {
    marker_gene_list <- read.csv(opt$gene_list, sep = ",", header=FALSE)
    marker_gene_list <- unlist(marker_gene_list["gene"], use.names = FALSE)
}

####################################################
### CALCULATE GINI SCORE AND SPEARMAN CORELATION ###
####################################################

# If no marker list is provided, run on ALL genes
eval_genes <- rownames(xenium_obj)

if (!is.null(marker_gene_list)){
    eval_genes <- intersect(rownames(xenium_obj), marker_gene_list)
}

# Create output dataframe
gene_pair_stats_df <- data.frame(
    sample = character(),
    gene1 = character(),
    gene2 = character(),
    spearman_r = numeric(),
    spearman_p = numeric(),
    gini_gene1 = numeric(),
    gini_gene2 = numeric(),
    filter_pass = logical(),
    stringsAsFactors = FALSE
)

for (i in seq_along(eval_genes)){

    for (j in seq_along(eval_genes)){

        if (i < j){ # Avoid duplicates (gene1-gene2 is the same as gene2-gene1)
            gene1 <- eval_genes[i]
            gene2 <- eval_genes[j]

            expr_data <- FetchData(xenium_obj, vars = c(gene1, gene2))

            # Skip evaluation if there's no expression for either gene
            if (all(expr_data[[gene1]] == 0) || all(expr_data[[gene2]] == 0)){
                next
            }

            # Spearman Correlation
            cor_test <- cor.test(
                expr_data[[gene1]],
                expr_data[[gene2]],
                method = "spearman"
            )

            # Gini Indices
            gini_gene1 <- calculate_gini_score(expr_data[[gene1]])
            gini_gene2 <- calculate_gini_score(expr_data[[gene2]])

            # Append result to dataframe
            gene_pair_stats_df <- gene_pair_stats_df %>% add_row(
                sample = unique(xenium_obj$Sample),
                gene1 = gene1,
                gene2 = gene2,
                spearman_r = cor_test$estimate,
                spearman_p = cor_test$p.value,
                gini_gene1 = gini_gene1,
                gini_gene2 = gini_gene2
            )
        }
    }
}

#################
### SAVE DATA ###
#################

write.csv(
    gene_pair_stats_df,
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
