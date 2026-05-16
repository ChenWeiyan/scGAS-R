#' Compute Single-Cell Gene Activation Scores via Network Propagation
#'
#' Propagates Metacell-level GAS to single-cell resolution in two steps:
#' \enumerate{
#'   \item \strong{Initialisation (scGAS0).}  Each non-centroid cell inherits a
#'     Gaussian-kernel weighted average of its \code{n_metacell_refs} nearest
#'     Metacell centroids in LSI space.
#'   \item \strong{Network propagation (scGAS).}  A mutual-KNN graph is built;
#'     the top \code{seed_fraction} cells seed a random-walk propagation.
#'     Scores are outlier-clipped and min-max scaled to \eqn{[0, 1]}.
#' }
#'
#' The resulting genes x cells scGAS matrix is stored as a new \code{"scGAS"}
#' assay in the returned \code{SeuratObject}.  PCA and UMAP on the scGAS
#' values are run when \code{run_dim_reduction = TRUE}.
#'
#' @param seurat_obj \code{SeuratObject} from \code{\link{scgas_preprocess}}.
#'   Must contain an \code{lsi} reduction, \code{mc_membership} metadata
#'   (from \code{\link{scgas_metacell}}), and fitted models in
#'   \code{seurat_obj@misc$scgas_models} (from
#'   \code{\link{scgas_train_models}}).
#' @param lsi_dims LSI dimensions to use.  Must match those used in
#'   \code{\link{scgas_preprocess}}.  Default \code{2:30}.
#' @param knn_k Mutual nearest neighbour count for graph construction.
#'   Default \code{30}.
#' @param n_metacell_refs Number of nearest Metacell centroids used to
#'   initialise each cell's scGAS0.  Default \code{3}.
#' @param seed_fraction Fraction of top-scoring cells used as propagation
#'   seeds.  Default \code{0.10}.
#' @param gamma Random-walk restart probability.  Default \code{0.05}.
#' @param outlier_quantile Upper quantile for outlier capping before scaling.
#'   Default \code{0.95}.
#' @param chunk_size Number of genes processed per parallel batch.  Default
#'   \code{500}.
#' @param run_dim_reduction Logical.  Run PCA and UMAP on the scGAS assay.
#'   Default \code{TRUE}.
#' @param n_pcs Number of PCs for PCA.  Default \code{30}.
#' @param n_neighbors \code{n.neighbors} for UMAP.  Default \code{30}.
#' @param n_cores Cores for \code{parallel::mclapply}.  Default \code{1}.
#' @param verbose Logical.  Default \code{TRUE}.
#' @param out_dir Base output directory.  \code{NULL} (default) inherits from
#'   \code{seurat_obj@@misc$out_dir}.  Pass an explicit path to override.
#' @param run_name Run label.  \code{NULL} (default) inherits from
#'   \code{seurat_obj@@misc$run_name}.  Pass an explicit string to override.
#'   The computed scGAS matrix is cached as
#'   \code{out_dir/run_name/scgas_mat.rds}; re-running with the same values
#'   loads the cached matrix and skips the computation.
#' @param save_obj Logical.  Save the returned \code{SeuratObject} to
#'   \code{out_dir/run_name/seurat_scgas_computed.rds}.  Default \code{FALSE}.
#'
#' @return The input \code{SeuratObject} with:
#'   \describe{
#'     \item{\code{scGAS} assay}{Genes x cells matrix of scGAS values in
#'       \eqn{[0, 1]}.}
#'     \item{\code{scgaspca}}{(If \code{run_dim_reduction}) PCA on scGAS.}
#'     \item{\code{scgasumap}}{(If \code{run_dim_reduction}) 2-D UMAP on
#'       scGAS PCA.}
#'   }
#'
#' @seealso \code{\link{scgas_add_assay}} to add a pre-computed scGAS matrix
#'   to a \code{SeuratObject} without re-running the computation.
#'
#' @examples
#' \dontrun{
#' obj <- scgas_compute(
#'   seurat_obj = obj,
#'   lsi_dims   = 2:30,
#'   knn_k      = 30,
#'   n_cores    = 8
#' )
#'
#' ## Visualise marker gene activation scores
#' FeaturePlot(obj, features = c("SLC17A7", "GAD1", "MBP"),
#'             reduction = "scgasumap", ncol = 3)
#' }
#'
#' @importFrom MatrixGenerics colMeans2
#' @importFrom parallel mclapply
#' @importFrom SCAVENGE getmutualknn randomWalk_sparse capOutlierQuantile
#'   max_min_scale
#' @importFrom Seurat CreateAssayObject DefaultAssay FindVariableFeatures
#'   ScaleData RunPCA RunUMAP VariableFeatures
#' @export
scgas_compute <- function(seurat_obj,
                          lsi_dims          = 2:30,
                          knn_k             = 30L,
                          n_metacell_refs   = 3L,
                          seed_fraction     = 0.10,
                          gamma             = 0.05,
                          outlier_quantile  = 0.95,
                          chunk_size        = 500L,
                          run_dim_reduction = TRUE,
                          n_pcs             = 30L,
                          n_neighbors       = 30L,
                          n_cores           = 1L,
                          verbose           = TRUE,
                          out_dir           = NULL,
                          run_name          = NULL,
                          save_obj          = FALSE) {

  stopifnot(inherits(seurat_obj, "Seurat"))

  ## ── Resolve run directory ─────────────────────────────────────────────────
  rd               <- .resolve_run_dir(out_dir, run_name, seurat_obj@misc, verbose)
  scgas_cache_path <- if (rd$use_cache) file.path(rd$run_dir, "scgas_mat.rds") else NULL

  if (rd$use_cache && file.exists(scgas_cache_path)) {
    if (verbose) message("[scGAS] Loading cached scGAS matrix from ", scgas_cache_path)
    scgas_mat <- readRDS(scgas_cache_path)
    seurat_obj@misc$out_dir  <- rd$out_dir
    seurat_obj@misc$run_name <- rd$run_name
    seurat_obj <- scgas_add_assay(seurat_obj, scgas_mat,
                                  run_dim_reduction, n_pcs, n_neighbors, verbose)
    if (save_obj) {
      obj_path <- file.path(rd$run_dir, "seurat_scgas_computed.rds")
      if (verbose) message("[scGAS] Saving computed object to ", obj_path)
      saveRDS(seurat_obj, obj_path)
    }
    return(seurat_obj)
  }

  if (!"lsi" %in% names(seurat_obj@reductions))
    stop("[scGAS] 'lsi' reduction not found. Run scgas_preprocess() first.")

  if (is.null(seurat_obj@misc$scgas_models))
    stop("[scGAS] No models found in seurat_obj@misc$scgas_models. ",
         "Run scgas_train_models() first.")

  if (!"mc_membership" %in% colnames(seurat_obj@meta.data))
    stop("[scGAS] 'mc_membership' not found in metadata. ",
         "Run scgas_metacell() first.")

  stopifnot(
    is.numeric(seed_fraction), seed_fraction > 0, seed_fraction <= 1,
    is.numeric(gamma), gamma > 0
  )

  trained_models <- seurat_obj@misc$scgas_models
  cl_labels      <- seurat_obj$mc_membership
  cl_levels      <- attr(trained_models, "cl_levels")
  sel_genes      <- names(trained_models)

  embeddings <- seurat_obj@reductions$lsi@cell.embeddings[, lsi_dims, drop = FALSE]
  n_cells    <- nrow(embeddings)

  if (verbose) message("[scGAS] Building mutual-KNN graph (k = ", knn_k, ")")
  mutualknn <- SCAVENGE::getmutualknn(embeddings, knn_k)

  if (verbose) message("[scGAS] Identifying Metacell centroids")
  centroid_cells <- .find_centroids(embeddings, cl_labels, cl_levels)

  if (verbose) message("[scGAS] Computing cell-to-centroid distances")
  centroid_emb      <- embeddings[centroid_cells, , drop = FALSE]
  dist_to_centroids <- .cross_dist(embeddings, centroid_emb)
  rownames(dist_to_centroids) <- rownames(embeddings)
  colnames(dist_to_centroids) <- centroid_cells

  k_band    <- min(knn_k, ncol(dist_to_centroids))
  alpha_vec <- apply(dist_to_centroids, 1, function(x) sort(x, partial = k_band)[k_band])
  names(alpha_vec) <- rownames(embeddings)

  other_cells <- base::setdiff(rownames(embeddings), centroid_cells)
  cells_ref   <- .build_cell_refs(other_cells, centroid_cells,
                                  dist_to_centroids, alpha_vec, n_metacell_refs)

  if (verbose) message("[scGAS] Computing scGAS for ", length(sel_genes), " genes")
  n_genes   <- length(sel_genes)
  n_chunks  <- max(1L, round(n_genes / chunk_size))
  gene_chunks <- split(seq_len(n_genes),
                       sample(factor(seq_len(n_genes) %% n_chunks)))

  results <- vector("list", n_chunks)
  for (chunk_i in seq_len(n_chunks)) {
    chunk_genes <- sel_genes[gene_chunks[[chunk_i]]]
    chunk_res   <- parallel::mclapply(
      X                = as.list(chunk_genes),
      FUN              = .compute_one_gene,
      trained_models   = trained_models,
      centroid_cells   = centroid_cells,
      other_cells      = other_cells,
      cells_ref        = cells_ref,
      n_cells          = n_cells,
      cell_names       = colnames(seurat_obj),
      mutualknn        = mutualknn,
      seed_fraction    = seed_fraction,
      gamma            = gamma,
      outlier_quantile = outlier_quantile,
      mc.cores         = n_cores
    )
    names(chunk_res)   <- chunk_genes
    results[[chunk_i]] <- chunk_res
    gc()
  }

  if (verbose) message("[scGAS] Assembling result matrix")
  scgas_mat <- matrix(
    NA_real_,
    nrow     = n_genes,
    ncol     = ncol(seurat_obj),
    dimnames = list(sel_genes, colnames(seurat_obj))
  )
  for (chunk_i in seq_len(n_chunks)) {
    for (g in names(results[[chunk_i]])) {
      res_g <- results[[chunk_i]][[g]]
      if (!is.null(res_g)) scgas_mat[g, ] <- res_g$scGAS
    }
  }

  if (verbose) message("[scGAS] Done: ", nrow(scgas_mat), " genes x ",
                       ncol(scgas_mat), " cells.")

  if (rd$use_cache) {
    dir.create(rd$run_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) message("[scGAS] Saving scGAS matrix to ", scgas_cache_path)
    saveRDS(scgas_mat, scgas_cache_path)
  }

  ## ── Embed scGAS matrix as an assay ───────────────────────────────────────
  seurat_obj <- scgas_add_assay(
    seurat_obj        = seurat_obj,
    scgas_mat         = scgas_mat,
    run_dim_reduction = run_dim_reduction,
    n_pcs             = n_pcs,
    n_neighbors       = n_neighbors,
    verbose           = verbose
  )

  seurat_obj@misc$out_dir  <- rd$out_dir
  seurat_obj@misc$run_name <- rd$run_name

  if (save_obj && rd$use_cache) {
    obj_path <- file.path(rd$run_dir, "seurat_scgas_computed.rds")
    if (verbose) message("[scGAS] Saving computed object to ", obj_path)
    saveRDS(seurat_obj, obj_path)
  }

  return(seurat_obj)
}


#' Add a Pre-Computed scGAS Matrix as an Assay Inside a SeuratObject
#'
#' Stores the genes x cells scGAS matrix as a new \code{"scGAS"} assay and
#' optionally runs PCA and UMAP on the scGAS values.  This function is
#' called internally by \code{\link{scgas_compute}}; use it directly only when
#' you have a pre-computed matrix and want to (re-)add it to the object.
#'
#' @param seurat_obj \code{SeuratObject}.
#' @param scgas_mat Genes x cells numeric matrix of scGAS values.
#' @param run_dim_reduction Logical.  Run PCA and UMAP on the scGAS assay.
#'   Default \code{TRUE}.
#' @param n_pcs Number of PCs for PCA and UMAP.  Default \code{30}.
#' @param n_neighbors \code{n.neighbors} for UMAP.  Default \code{30}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return The \code{SeuratObject} with a new \code{scGAS} assay and
#'   optionally \code{scgaspca} and \code{scgasumap} reductions.
#'
#' @examples
#' \dontrun{
#' ## Re-run dim reduction with different parameters without recomputing GAS
#' obj <- scgas_add_assay(obj, scgas_mat = as.matrix(obj@assays$scGAS@data),
#'                        n_pcs = 20, n_neighbors = 15)
#' }
#'
#' @importFrom Seurat CreateAssayObject DefaultAssay FindVariableFeatures
#'   ScaleData RunPCA RunUMAP VariableFeatures
#' @export
scgas_add_assay <- function(seurat_obj,
                            scgas_mat,
                            run_dim_reduction = TRUE,
                            n_pcs             = 30L,
                            n_neighbors       = 30L,
                            verbose           = TRUE) {
  stopifnot(inherits(seurat_obj, "Seurat"), is.matrix(scgas_mat))

  assay_obj      <- Seurat::CreateAssayObject(counts = scgas_mat)
  assay_obj@data <- scgas_mat
  seurat_obj[["scGAS"]] <- assay_obj

  if (run_dim_reduction) {
    if (verbose) message("[scGAS] Running dimensionality reduction on scGAS assay")
    Seurat::DefaultAssay(seurat_obj) <- "scGAS"
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = verbose)
    seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = verbose)
    seurat_obj <- Seurat::RunPCA(
      seurat_obj,
      features       = Seurat::VariableFeatures(seurat_obj),
      reduction.name = "scgaspca",
      verbose        = verbose
    )
    seurat_obj <- Seurat::RunUMAP(
      seurat_obj,
      reduction      = "scgaspca",
      dims           = seq_len(n_pcs),
      n.components   = 2L,
      n.neighbors    = n_neighbors,
      reduction.name = "scgasumap",
      verbose        = verbose
    )
  }
  return(seurat_obj)
}


## ─────────────────────────────────────────────────────────────────────────────
## Internal helpers
## ─────────────────────────────────────────────────────────────────────────────

.cross_dist <- function(A, B) {
  sq_a <- rowSums(A^2)
  sq_b <- rowSums(B^2)
  cp   <- tcrossprod(A, B)
  sqrt(pmax(outer(sq_a, sq_b, "+") - 2 * cp, 0))
}

.find_centroids <- function(embeddings, cl_labels, cl_levels) {
  centroids <- character(length(cl_levels))
  names(centroids) <- paste0("Cluster-", cl_levels)
  for (cl in cl_levels) {
    sel_cells <- which(cl_labels == cl)
    cl_emb    <- embeddings[sel_cells, , drop = FALSE]
    cl_mean   <- MatrixGenerics::colMeans2(cl_emb)
    dists     <- apply(cl_emb, 1, function(x) sqrt(sum((x - cl_mean)^2)))
    centroids[paste0("Cluster-", cl)] <- names(sort(dists))[1]
  }
  return(centroids)
}

.build_cell_refs <- function(other_cells, centroid_cells,
                             dist_to_centroids, alpha_vec, n_metacell_refs) {
  refs <- vector("list", length(other_cells))
  names(refs) <- other_cells
  for (cell in other_cells) {
    alpha     <- alpha_vec[cell]
    raw_dists <- dist_to_centroids[cell, centroid_cells]
    weights   <- exp(-(raw_dists / alpha)^2)
    weights   <- weights / sum(weights)
    top_idx   <- order(weights, decreasing = TRUE)[seq_len(min(n_metacell_refs, length(weights)))]
    top_wgt   <- weights[top_idx]
    refs[[cell]] <- top_wgt / sum(top_wgt)
  }
  return(refs)
}

.compute_one_gene <- function(gene_name,
                              trained_models,
                              centroid_cells,
                              other_cells,
                              cells_ref,
                              n_cells,
                              cell_names,
                              mutualknn,
                              seed_fraction,
                              gamma,
                              outlier_quantile) {
  tryCatch({
    model_g <- trained_models[[gene_name]]
    gas_mc  <- model_g$GAS
    names(gas_mc) <- unname(centroid_cells)

    scgas0 <- stats::setNames(numeric(length(cell_names)), cell_names)
    scgas0[unname(centroid_cells)] <- gas_mc

    for (cell in other_cells) {
      ref_wgt      <- cells_ref[[cell]]
      scgas0[cell] <- sum(ref_wgt * gas_mc[names(ref_wgt)])
    }

    seed_n   <- max(1L, floor(seed_fraction * length(scgas0)))
    seed_idx <- rownames(mutualknn)[order(scgas0, decreasing = TRUE)[seq_len(seed_n)]]
    np_score <- SCAVENGE::randomWalk_sparse(intM = mutualknn,
                                            queryCells = seed_idx,
                                            gamma = gamma)
    np_score2 <- SCAVENGE::capOutlierQuantile(np_score, outlier_quantile)
    scgas     <- as.numeric(SCAVENGE::max_min_scale(np_score2))
    names(scgas) <- rownames(mutualknn)
    scgas <- scgas[cell_names]

    list(scGAS0 = scgas0, scGAS = scgas)
  }, error = function(e) NULL)
}
