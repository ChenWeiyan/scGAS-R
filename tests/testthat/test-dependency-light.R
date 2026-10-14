# Tests for the parts of the pipeline that do not need fragment files or the
# ENCODE reference (metacell construction, the compute-step linear-algebra
# helpers, assay attachment, the Chromatin Potential Field, and cpf S3
# dispatch). The heavy fragment -> train -> compute path is exercised by
# toy_test/run_toy_pipeline.R at the project level.

suppressWarnings(suppressMessages({
  library(Seurat)
  library(Matrix)
  library(igraph)
}))

# ---- helpers ---------------------------------------------------------------

make_block_graph_obj <- function(nclust = 4, per = 25, seed = 1) {
  set.seed(seed)
  n     <- nclust * per
  cells <- paste0("cell", seq_len(n))
  blk   <- kronecker(diag(nclust), matrix(1, per, per))
  noise <- matrix(runif(n * n) < 0.02, n, n)
  A <- (blk + noise) > 0
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  A <- Matrix(A * 1.0, sparse = TRUE)
  dimnames(A) <- list(cells, cells)

  cnt <- Matrix(rpois(30 * n, 2), 30, n, sparse = TRUE)
  rownames(cnt) <- paste0("g", seq_len(30))
  colnames(cnt) <- cells
  obj <- CreateSeuratObject(counts = cnt, assay = "ATAC")
  obj@graphs[["ATAC_snn"]] <- as.Graph(A)
  obj
}

# ---- scgas_metacell --------------------------------------------------------

test_that("scgas_metacell partitions a block graph into metacells", {
  obj <- make_block_graph_obj()
  obj <- scgas_metacell(obj, gamma = 12, min_size = 5, max_size = 40,
                        verbose = FALSE)

  expect_true("mc_membership" %in% colnames(obj@meta.data))
  expect_length(obj$mc_membership, ncol(obj))
  expect_false(is.null(obj@misc$mc_membership))
  # labels are a contiguous 1:K relabelling
  expect_identical(sort(unique(obj$mc_membership)),
                   seq_len(max(obj$mc_membership)))
})

# ---- compute-step helpers --------------------------------------------------

test_that("compute helpers return well-formed output", {
  set.seed(2)
  n     <- 100
  cells <- paste0("cell", seq_len(n))
  emb   <- matrix(rnorm(n * 5), n, 5,
                  dimnames = list(cells, paste0("LSI_", 1:5)))
  cl_labels <- rep(1:4, each = 25)
  cl_levels <- levels(as.factor(cl_labels))

  cent <- scGAS:::.find_centroids(emb, cl_labels, cl_levels)
  expect_length(cent, length(cl_levels))
  expect_true(all(cent %in% cells))

  d <- scGAS:::.cross_dist(emb, emb[cent, , drop = FALSE])
  expect_equal(dim(d), c(n, length(cent)))
  expect_true(all(d >= 0))

  alpha <- apply(d, 1, function(x) sort(x, partial = 2)[2])
  names(alpha) <- cells
  rownames(d) <- cells; colnames(d) <- cent
  others <- setdiff(cells, unname(cent))
  W <- scGAS:::.build_weight_matrix(others, unname(cent), d, alpha, 3)
  expect_equal(nrow(W), length(others))
  # each cell's weights over its reference centroids sum to 1
  expect_true(all(abs(Matrix::rowSums(W) - 1) < 1e-8))

  A  <- Matrix(abs(rsparsematrix(n, n, 0.1)), sparse = TRUE)
  M  <- A / pmax(Matrix::colSums(A), 1)
  P0 <- matrix(0, n, 2); P0[1:5, 1] <- 0.2; P0[6:10, 2] <- 0.2
  P  <- scGAS:::.batched_rw(M, P0, gamma = 0.05)
  expect_equal(dim(P), dim(P0))
  expect_true(all(is.finite(P)))
})

# ---- scgas_add_assay -------------------------------------------------------

test_that("scgas_add_assay attaches a scGAS assay", {
  obj   <- make_block_graph_obj()
  genes <- rownames(obj)
  scg   <- matrix(runif(length(genes) * ncol(obj)), length(genes), ncol(obj),
                  dimnames = list(genes, colnames(obj)))
  obj   <- scgas_add_assay(obj, scg, run_dim_reduction = FALSE, verbose = FALSE)
  expect_true("scGAS" %in% names(obj@assays))
})

# ---- scgas_chromatin_potential + cpf S3 methods ----------------------------

make_cpf <- function() {
  obj   <- make_block_graph_obj()
  genes <- rownames(obj)
  cells <- colnames(obj)
  scg   <- matrix(runif(length(genes) * length(cells)),
                  length(genes), length(cells), dimnames = list(genes, cells))
  obj   <- scgas_add_assay(obj, scg, run_dim_reduction = FALSE, verbose = FALSE)
  obj[["RNA"]] <- CreateAssayObject(
    counts = matrix(abs(rnorm(length(genes) * length(cells))),
                    length(genes), length(cells),
                    dimnames = list(genes, cells)))
  obj@reductions[["scgasumap"]] <- CreateDimReducObject(
    embeddings = matrix(rnorm(length(cells) * 2), length(cells), 2,
                        dimnames = list(cells, c("UMAP_1", "UMAP_2"))),
    key = "scgasumap_", assay = "scGAS")
  scgas_chromatin_potential(obj, n_hvg = 25, knn_k = 5, verbose = FALSE)
}

test_that("scgas_chromatin_potential returns a well-formed cpf object", {
  cpf <- make_cpf()
  expect_s3_class(cpf, "cpf")
  expect_true(all(c("cpf_df", "plot", "hvg", "scgas_hvg", "rna_hvg",
                    "embedding", "cell_metadata", "params") %in% names(cpf)))
  expect_s3_class(cpf$plot, "ggplot")
  expect_equal(nrow(cpf$cpf_df), 100)
})

test_that("cpf S3 methods dispatch (guards NAMESPACE S3method registration)", {
  cpf <- make_cpf()
  # If print.cpf were merely export()ed instead of registered via
  # S3method(print, cpf), these would fall through to print.default.
  expect_output(print(cpf), "Chromatin Potential Field")
  expect_output(summary(cpf), "Arrow magnitude")
})
