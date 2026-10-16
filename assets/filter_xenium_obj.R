
library(Seurat)     # Main analysis package

#' This function filters a Seurat object processed by the nf_xpatial pipeline
#' 
#' @param input The Seurat object to filter
#' @param assays A character vector of assays to keep. Values: ("Xenium", "AreaNorm")
#' @param clustering_method A character vector of clustering methods to keep. Values ("Seurat", "BANKSY")
#' @param srt_dims A double vector of dimensions used for Seurat clustering to keep. Leave empty to include all dimensions.
#' @param srt_res A double vector of resolutions used for Seurat clustering to keep. Leave empty to include all resolutions.
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
    srt_dims = vector(),
    srt_res = vector(),
    bsky_dims = vector(),
    bsky_res = vector(),
    bsky_lambda = vector(),
    bsky_kgeom = vector() ) {

  metadata_cols <- colnames(input@meta.data)
  cluster_cols_keep <- metadata_cols[!grepl("^clust_", metadata_cols)]
  
  reductions_keep <- c()
  
  # ------------------------ #
  # --- HELPER FUNCTIONS --- #
  # ------------------------ #
  
  # Append "<prefix><value>" for every value (or ".*" if no values supplied)
  add_param_to_keep <- function(prefix, param_list = vector() ) {
    suffix <- ifelse(length(param_list), param_list, ".*" )
    return(
      apply(expand.grid(prefix, suffix), 1, function(x) paste(x, collapse = ""))
    )
  }
  
  # Seurat-style clustering: _d<dims>_r<res>
  seurat_patterns <- function(clust_prefix, reduction_prefix) {
    clusts <- clust_prefix |>
      paste0("_d") |> add_param_to_keep(srt_dims) |>
      paste0("_r") |> add_param_to_keep(srt_res)

    return (
      list(
        clusters   = clusts,
        reductions = add_param_to_keep(paste0(reduction_prefix, "_.*_d"), srt_dims)
      )
    )
  }

  # BANKSY-style clustering: _l<lambda>_k<k_geom>_d<dims>_r<res>
  banksy_patterns <- function(clust_prefix, reduction_prefix) {
    clusts <- clust_prefix |>
      paste0("_l") |> add_param_to_keep(bsky_lambda) |>
      paste0("_k") |> add_param_to_keep(bsky_kgeom) |>
      paste0("_d") |> add_param_to_keep(bsky_dims) |>
      paste0("_r") |> add_param_to_keep(bsky_res)

    list(
      clusters   = clusts,
      reductions = add_param_to_keep(paste0(reduction_prefix, "_.*_d"), bsky_dims)
    )
  }
  
  # Define the patterns
  method_patterns <- list(
    Seurat       = function() seurat_patterns("clust_SEU",     "Seurat"),
    BANKSYSeurat = function() banksy_patterns("clust_BSKYSEU", "BANKSYSeurat"),
    BANKSY       = function() banksy_patterns("clust_BSKY",    "BANKSY")
  )
  
  # ---------------------------- #
  # --- COLLECT WHAT TO KEEP --- #
  # ---------------------------- #
  
  cluster_cols_keep <- metadata_cols[!grepl("^clust_", metadata_cols)]
  reductions_keep   <- c()
  
  for (method in intersect(names(method_patterns), clustering_method)) {
    pats <- method_patterns[[method]]()
    cluster_cols_keep <- c(cluster_cols_keep, pats$clusters)
    reductions_keep   <- c(reductions_keep, pats$reductions)
  }
  
  cluster_keep_regex    <- paste(cluster_cols_keep, collapse = "|")
  reductions_keep_regex <- paste(reductions_keep, collapse = "|")
  
  
  # -------------- #
  # --- FILTER --- #
  # -------------- #
  filtered_seurat <- input
  
  all_reductions <- Reductions(input)
  reductions_to_remove <- if (length(reductions_keep)) {
    all_reductions[!grepl(reductions_keep_regex, all_reductions)]
  } else {
    all_reductions
  }
  for (reduction in reductions_to_remove) {
    filtered_seurat[[reduction]] <- NULL
  }
  
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
#  bsky_lambda = c(0.9),
#  bsky_kgeom = c(15),
#  bsky_dims = c(30),
#  bsky_res = c(0.8)
#)

# This will grab only the Seurat cluster and reductions for dimension 20 and resolution 0.4
#single_seurat <- filter_xenium_obj(
#  input = full_object,
#  clustering_method = c("Seurat"),
#  srt_dims = c(20),
#  srt_res = c(0.4)
#)

# This will grab only the Seurat cluster and reductions for dimension 20, 25,and 30 and resolution 0.4
#mutliple_dims_seurat <- filter_xenium_obj(
#   input = full_object,
#   clustering_method = c("Seurat"),
#   srt_dims = c(20,25,30),
#   srt_res = c(0.4)
# )

# This will grab BanksySeurat cluster and reductions for lambda 0.9, k_geom 15, nPCs 30, and resolution 0.8
#single_banksy <- filter_xenium_obj(
#  input = full_object,
#  clustering_method = c("BANKSYSeurat"),
#  bsky_lambda = c(0.9),
#  bsky_kgeom = c(15),
#  bsky_dims = c(30),
#  bsky_res = c(0.8)
#)