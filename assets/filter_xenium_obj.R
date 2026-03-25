
library(Seurat)     # Main analysis package

#' This function filters a Seurat object processed by the nf_xpatial pipeline
#' 
#' @param input The Seurat object to filter
#' @param assays A character vector of assays to keep. Values: ("Xenium", "AreaNorm")
#' @param clustering_method A character vector of clustering methods to keep. Values ("Seurat", "BANKSY")
#' @param hmy_dims A double vector of dimensions used for Seurat clustering to keep. Leave empty to include all dimensions.
#' @param hmy_res A double vector of resolutions used for Seurat clustering to keep. Leave empty to include all resolutions.
#' @param bsky_dims A double vector of dimensions used for BANKSY clustering to keep. Leave empty to include all dimensions.
#' @param bsky_res A double vector of resolutions used for BANKSY clustering to keep. Leave empty to include all resolutions.
#' @param bsky_lambda A double vector of lambdas used for BANKSY clustering to keep. Leave empty to include all lambdas.
#' @param bsky_kgeom A double vector of kgeom used for BANKSY clustering to keep. Leave empty to include all k_geoms.
#'
#' @returns A Seurat object containing only the assays, reductions, and clusters listed above.

filter_xenium_obj <- function(
    input = NULL,
    assays = "Xenium",
    clustering_method = vector(),
    hmy_dims = vector(),
    hmy_res = vector(),
    bsky_dims = vector(),
    bsky_res = vector(),
    bsky_lambda = vector(),
    bsky_kgeom = vector() ) {

    metadata_cols <- colnames(input@meta.data)
    cluster_cols_keep <- metadata_cols[!grepl("^clust_", metadata_cols)]

    reductions_keep <- c()

    # Helper: build parameter patterns for regex matching
    add_param_to_keep <- function(prefix, param_list = vector() ) {

        suffix <- ifelse(length(param_list), param_list, ".*" )

        return(
          apply(expand.grid(prefix, suffix), 1, function(x) paste(x, collapse = ""))
        )
    }

    if ("Seurat" %in% clustering_method) {
      # Get the Cluster columns
      srt_clusts <- c("clust_HMY")

      srt_clusts <- add_param_to_keep(paste0(srt_clusts, "_d"), hmy_dims)
      srt_clusts <- add_param_to_keep(paste0(srt_clusts, "_r"), hmy_res)

      cluster_cols_keep <- c(cluster_cols_keep, srt_clusts)

      # Get the reductions
      reductions_keep <- c(
        reductions_keep,
        add_param_to_keep("Harmony_.*_d", hmy_dims)
      )
    }

    if ("BANKSY" %in% clustering_method) {
      # Get the Cluster columns
      bsky_clusts <- c("clust_BSKY_AGF1")

      bsky_clusts <- add_param_to_keep(paste0(bsky_clusts, "_L"), bsky_lambda)
      bsky_clusts <- add_param_to_keep(paste0(bsky_clusts, "_k"), bsky_kgeom)
      bsky_clusts <- add_param_to_keep(paste0(bsky_clusts, "_d"), bsky_dims)
      bsky_clusts <- add_param_to_keep(paste0(bsky_clusts, "_R"), bsky_res)

      cluster_cols_keep <- c(cluster_cols_keep, bsky_clusts)

      # Get the reductions
      reductions_keep <- c(
        reductions_keep,
        add_param_to_keep("BANKSY_.*_d", bsky_dims)
      )
    }

    filtered_seurat <- input

    # Build regex patterns
    reductions_keep_regex <- paste(reductions_keep, collapse = "|")
    cluster_keep_regex <- paste(cluster_cols_keep, collapse = "|")

    # Filter reductions
    all_reductions <- Reductions(input)
    reductions_to_remove <- all_reductions[!grepl(reductions_keep_regex, all_reductions)]

    for (reduction in reductions_to_remove) {
      filtered_seurat[[reduction]] <- NULL
    }

    # Filter metadata columns
    cols_to_remove <- metadata_cols[!grepl(cluster_keep_regex, metadata_cols)]

    for (col in cols_to_remove) {
      filtered_seurat[[col]] <- NULL
    }

    return(filtered_seurat)

  }

### Examples ###
# This will grab all the BANKSY cluster and reductions, filtering out the Seurat cluster results
#only_banksy_obj <- filter_xenium_obj(
#  input = full_object,
#  clustering_method = c("BANKSY")
#)

# This will grab the BANKSY cluster and reductions for lambda 0.9, k_geom 15, nPCs 30, and resolution 0.8
#single_banksy <- filter_xenium_obj(
#  input = full_object,
#  clustering_method = c("Seurat"),
#  bksy_lambda = c(0.9),
#  bsky_kgeom = c(15),
#  bsky_dims = c(30),
#  bsky_res = c(0.8)
#)

# This will grab only the Seurat cluster and reductions for dimension 20 and resolution 0.4
#single_seurat <- filter_xenium_obj(
#  input = full_object,
#  clustering_method = c("Seurat"),
#  hmy_dims = c(20),
#  hmy_res = c(0.4)
#)

# This will grab only the Seurat cluster and reductions for dimension 20, 25,and 30 and resolution 0.4
#mutliple_dims_seurat <- filter_xenium_obj(
#   input = full_object,
#   clustering_method = c("Seurat"),
#   hmy_dims = c(20,25,30),
#   hmy_res = c(0.4)
# )