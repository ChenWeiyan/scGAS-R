#' Compute Single-Cell Gene Activation Scores via Network Propagation
#'
#' Propagates Metacell-level GAS to single-cell resolution in two steps:
#' \enumerate{
#'   \item \strong{Initialisation (scGAS0).}  Each non-centroid cell inherits a
#'     Gaussian-kernel weighted average of its \code{n_metacell_refs} nearest
#'     Metacell centroids in LSI space.  The kernel bandwidth of each cell is
#'     its distance to its \code{kernel_k}-th nearest neighbour among
#'     \emph{all} cells, so the kernel adapts to local cell density.
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
#' @param kernel_k Integer.  Neighbour rank defining the per-cell Gaussian
#'   kernel bandwidth: \eqn{\alpha_i} is the LSI distance from cell \eqn{i}
#'   to its \code{kernel_k}-th nearest other cell (self excluded), and the
#'   weight to centroid \eqn{m} is \eqn{\exp(-(d_{im}/\alpha_i)^2)}.  Smaller
#'   values give a sharper, more local kernel.  Distances are computed exactly
#'   in blocks, so memory is \eqn{O(n)} per block rather than
#'   \eqn{O(n^2)}.  Default \code{30}.
#' @param seed_fraction Fraction of top-scoring cells used as propagation
#'   seeds.  Default \code{0.10}.
#' @param gamma Random-walk restart probability.  Default \code{0.05}.
#' @param outlier_quantile Upper quantile for outlier capping before scaling.
#'   Default \code{0.95}.
#' @param chunk_size Number of genes per memory batch.  All genes in a batch
#'   are scored simultaneously via batched sparse matrix algebra; larger values
#'   are faster but use more RAM.  Default \code{500}.
#' @param run_dim_reduction Logical.  Run PCA and UMAP on the scGAS assay.
#'   Default \code{TRUE}.
#' @param n_pcs Number of PCs for PCA.  Default \code{30}.
#' @param n_neighbors \code{n.neighbors} for UMAP.  Default \code{30}.
#' @param n_cores Kept for backward compatibility.  Parallelism is provided by
#'   the BLAS/OpenMP threads of the linked matrix library (set
#'   \code{OMP_NUM_THREADS} before starting R).  Default \code{1}.
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
#' @importFrom Matrix sparseMatrix
#' @importFrom MatrixGenerics colMeans2
#' @importFrom SCAVENGE getmutualknn
#' @importFrom Seurat CreateAssayObject DefaultAssay FindVariableFeatures
#'   ScaleData RunPCA RunUMAP VariableFeatures
#' @export
scgas_compute <- function(seurat_obj,
                          lsi_dims          = 2:30,
                          knn_k             = 30L,
                          n_metacell_refs   = 3L,
                          kernel_k          = 30L,
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
    is.numeric(kernel_k), length(kernel_k) == 1L, kernel_k >= 1,
    is.numeric(seed_fraction), seed_fraction > 0, seed_fraction <= 1,
    is.numeric(gamma), gamma > 0
  )

  trained_models <- seurat_obj@misc$scgas_models
  cl_labels      <- seurat_obj$mc_membership
  cl_levels      <- attr(trained_models, "cl_levels")
  sel_genes      <- names(trained_models)
  cell_names     <- colnames(seurat_obj)

  embeddings <- seurat_obj@reductions$lsi@cell.embeddings[, lsi_dims, drop = FALSE]

  if (verbose) message("[scGAS] Building mutual-KNN graph (k = ", knn_k, ")")
  mutualknn <- SCAVENGE::getmutualknn(embeddings, knn_k)
  ## Column-stochastic transition matrix, as in SCAVENGE::randomWalk_sparse
  ## (intM <- t(t(intM) / colSums(intM))).  Built once and reused for all genes.
  deg       <- Matrix::colSums(mutualknn)
  trans_mat <- mutualknn %*% Matrix::Diagonal(x = 1 / pmax(deg, 1))
  dimnames(trans_mat) <- dimnames(mutualknn)
  rm(mutualknn)

  if (verbose) message("[scGAS] Identifying Metacell centroids")
  centroid_cells <- .find_centroids(embeddings, cl_labels, cl_levels)

  if (verbose) message("[scGAS] Computing cell-to-centroid distances")
  centroid_emb      <- embeddings[centroid_cells, , drop = FALSE]
  dist_to_centroids <- .cross_dist(embeddings, centroid_emb)
  rownames(dist_to_centroids) <- rownames(embeddings)
  colnames(dist_to_centroids) <- centroid_cells

  ## Per-cell bandwidth: distance to the kernel_k-th nearest cell among ALL
  ## cells (local scale), not among centroids.
  kernel_k <- min(as.integer(kernel_k), nrow(embeddings) - 1L)
  if (verbose) message("[scGAS] Computing Gaussian kernel bandwidth (kernel_k = ",
                       kernel_k, ")")
  alpha_vec <- .knn_bandwidth(embeddings, kernel_k)

  centroid_names <- unname(centroid_cells)
  centroid_idx   <- match(centroid_names, cell_names)
  other_cells    <- base::setdiff(cell_names, centroid_names)
  other_idx      <- match(other_cells, cell_names)

  W_init <- .build_weight_matrix(other_cells, centroid_names,
                                 dist_to_centroids, alpha_vec, n_metacell_refs)
  rm(dist_to_centroids); gc()

  ## Stack all per-gene metacell GAS vectors into one matrix (n_centroids × n_genes)
  gas_mc_mat <- do.call(cbind, lapply(trained_models, function(m) m$GAS))
  colnames(gas_mc_mat) <- sel_genes

  n_genes     <- length(sel_genes)
  n_seed      <- max(1L, floor(seed_fraction * length(cell_names)))
  n_chunks    <- max(1L, round(n_genes / chunk_size))
  gene_chunks <- split(seq_len(n_genes),
                       sample(factor(seq_len(n_genes) %% n_chunks)))

  if (verbose)
    message("[scGAS] Computing scGAS for ", n_genes, " genes in ",
            n_chunks, " batch", if (n_chunks > 1) "es" else "",
            " of ~", chunk_size, " genes (batched sparse algebra)")

  scgas_mat  <- matrix(
    NA_real_,
    nrow     = n_genes,
    ncol     = length(cell_names),
    dimnames = list(sel_genes, cell_names)
  )
  genes_done <- 0L

  for (chunk_i in seq_len(n_chunks)) {
    chunk_idx   <- gene_chunks[[chunk_i]]
    chunk_genes <- sel_genes[chunk_idx]
    nc          <- length(chunk_genes)

    ## ── scgas0: metacell GAS → single-cell initialisation ──────────────────
    gmc <- gas_mc_mat[, chunk_idx, drop = FALSE]            # n_centroids × nc
    s0  <- matrix(0, length(cell_names), nc,
                  dimnames = list(cell_names, chunk_genes))
    s0[centroid_idx, ] <- gmc
    s0[other_idx,    ] <- as.matrix(W_init %*% gmc)        # single matmul

    ## ── Seed matrix P0: top seed_fraction cells per gene ───────────────────
    ## Exactly n_seed seeds per gene (ties broken by cell order, as in the
    ## original sort(..., decreasing = TRUE, index.return = TRUE)), each with
    ## uniform weight 1 / n_seed.  Note: sort(partial =, decreasing = TRUE) is
    ## not supported by R and must not be used here.
    seed_idx <- apply(s0, 2, function(x) order(x, decreasing = TRUE)[seq_len(n_seed)])
    seed_idx <- matrix(seed_idx, nrow = n_seed)             # keep matrix shape if nc == 1
    P0       <- matrix(0, length(cell_names), nc)
    P0[cbind(as.vector(seed_idx), rep(seq_len(nc), each = n_seed))] <- 1 / n_seed

    ## ── Batched random walk ─────────────────────────────────────────────────
    ## Replaces nc individual SCAVENGE::randomWalk_sparse calls with one
    ## sparse × dense matrix iteration: trans_mat %*% P  (BLAS3)
    P <- .batched_rw(trans_mat, P0, gamma)
    rm(gmc, s0, seed_idx, P0); gc()

    ## ── Outlier cap + min-max scale per gene (column of P) ─────────────────
    caps <- apply(P, 2, quantile, probs = outlier_quantile, names = FALSE)
    for (j in seq_len(nc)) P[, j] <- pmin(P[, j], caps[j])
    mins <- apply(P, 2, min,      na.rm = TRUE)
    maxs <- apply(P, 2, max,      na.rm = TRUE)
    rng  <- pmax(maxs - mins, .Machine$double.eps)
    P_sc <- sweep(sweep(P, 2, mins, "-"), 2, rng, "/")

    scgas_mat[chunk_genes, ] <- t(P_sc)
    rm(P, P_sc); gc()

    genes_done <- genes_done + nc
    if (verbose) .scgas_progress(genes_done, n_genes)
  }

  if (verbose) {
    cat("\n", file = stderr())
    message("[scGAS] Done: ", nrow(scgas_mat), " genes x ",
            ncol(scgas_mat), " cells.")
  }

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

  assay_obj <- Seurat::CreateAssayObject(counts = scgas_mat)
  if (inherits(assay_obj, "Assay5")) {
    assay_obj <- Seurat::SetAssayData(assay_obj, layer = "data", new.data = scgas_mat)
  } else {
    assay_obj@data <- scgas_mat
  }
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
    actual_pcs <- ncol(seurat_obj@reductions$scgaspca@cell.embeddings)
    umap_dims  <- seq_len(min(n_pcs, actual_pcs))
    seurat_obj <- Seurat::RunUMAP(
      seurat_obj,
      reduction      = "scgaspca",
      dims           = umap_dims,
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

.scgas_progress <- function(done, total, width = 30L) {
  pct    <- done / total
  filled <- round(width * pct)
  bar    <- paste0(strrep("=", filled), strrep(" ", width - filled))
  cat(sprintf("\r[scGAS] [%s] %d/%d genes (%.0f%%)",
              bar, done, total, 100 * pct),
      file = stderr())
}

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

## Per-cell Gaussian kernel bandwidth: Euclidean distance from each cell to
## its k-th nearest OTHER cell (self excluded).  Exact all-cell distances are
## computed block by block, so peak memory is block_size x n_cells instead of
## the full n_cells x n_cells matrix.
.knn_bandwidth <- function(embeddings, k, block_size = 1000L) {
  n <- nrow(embeddings)
  stopifnot(k >= 1L, k <= n - 1L)
  alpha  <- numeric(n)
  starts <- seq.int(1L, n, by = block_size)
  for (s in starts) {
    idx <- s:min(s + block_size - 1L, n)
    blk <- .cross_dist(embeddings[idx, , drop = FALSE], embeddings)
    ## position 1 after sorting is the cell itself (distance 0) -> take k + 1
    alpha[idx] <- apply(blk, 1, function(x) sort(x, partial = k + 1L)[k + 1L])
  }
  names(alpha) <- rownames(embeddings)
  ## guard against duplicated embeddings (zero bandwidth)
  pmax(alpha, .Machine$double.eps)
}

## Builds a sparse n_other x n_centroids weight matrix.
## Row i holds the normalised Gaussian weights for other_cells[i] over its
## top n_ref nearest Metacell centroids.  Replaces the per-cell list built by
## .build_cell_refs, allowing vectorised W_init %*% gas_mc in workers.
## The kernel is monotone in distance, so the top n_ref weights are the n_ref
## nearest centroids.  The exponent is shifted by its minimum before exp() so
## that a small bandwidth cannot underflow all weights to 0 (-> NaN); the
## shift cancels in the normalisation.
.build_weight_matrix <- function(other_cells, centroid_names,
                                 dist_to_centroids, alpha_vec, n_metacell_refs) {
  n_other <- length(other_cells)
  n_ref   <- min(n_metacell_refs, length(centroid_names))
  ri <- integer(n_other * n_ref)
  ci <- integer(n_other * n_ref)
  xv <- numeric(n_other * n_ref)

  for (i in seq_along(other_cells)) {
    cell      <- other_cells[i]
    alpha     <- alpha_vec[cell]
    raw_dists <- dist_to_centroids[cell, centroid_names]
    top_j     <- order(raw_dists)[seq_len(n_ref)]
    z         <- (raw_dists[top_j] / alpha)^2
    wgt_top   <- exp(-(z - min(z)))
    wgt_top   <- wgt_top / sum(wgt_top)
    idx       <- seq.int((i - 1L) * n_ref + 1L, length.out = n_ref)
    ri[idx]   <- i
    ci[idx]   <- top_j
    xv[idx]   <- wgt_top
  }

  Matrix::sparseMatrix(
    i        = ri, j = ci, x = xv,
    dims     = c(n_other, length(centroid_names)),
    dimnames = list(other_cells, centroid_names)
  )
}

## Batched random walk with restart.
## M    — column-stochastic sparse transition matrix (n_cells × n_cells)
## P0   — dense seed matrix (n_cells × n_genes), columns sum to 1
## Returns a dense matrix of the same dimensions as P0.
## Convergence: iterate until the max per-gene L1 change < tol (SCAVENGE rule), or max_iter.
.batched_rw <- function(M, P0, gamma, max_iter = 1000L, tol = 1e-5) {
  c1 <- 1 - gamma
  P  <- P0                                 # base::matrix throughout
  for (i in seq_len(max_iter)) {
    P_new <- c1 * as.matrix(M %*% P) + gamma * P0   # coerce dgeMatrix → matrix
    ## same stopping rule as SCAVENGE::randomWalk_sparse (L1 change per gene)
    delta <- max(colSums(abs(P_new - P)))
    P <- P_new
    if (delta < tol) break
  }
  if (delta >= tol)
    warning("[scGAS] Random walk did not converge in ", max_iter,
            " iterations (max L1 change = ", signif(delta, 3), ").")
  P
}
