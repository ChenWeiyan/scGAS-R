## ─────────────────────────────────────────────────────────────────────────────
## cpf S3 class
## ─────────────────────────────────────────────────────────────────────────────

#' Chromatin Potential Field object
#'
#' An S3 class that stores all data needed to inspect and re-plot the Chromatin
#' Potential Field, independent of any \code{SeuratObject}.
#'
#' @section Slots:
#' \describe{
#'   \item{\code{cpf_df}}{Data frame with one row per cell: \code{x}, \code{y}
#'     (embedding position), \code{dx}, \code{dy} (arrow displacement),
#'     \code{magnitude} (arrow length), \code{cell} (barcode), and any
#'     metadata column requested via \code{point_colour}.}
#'   \item{\code{plot}}{\code{ggplot2} object showing background cells and CPF
#'     arrows.}
#'   \item{\code{hvg}}{Character vector of highly variable genes used for the
#'     Spearman correlation.}
#'   \item{\code{scgas_hvg}}{genes x cells scGAS matrix restricted to
#'     \code{hvg}.}
#'   \item{\code{rna_hvg}}{genes x cells log-normalised RNA matrix restricted
#'     to \code{hvg}.}
#'   \item{\code{embedding}}{cells x 2 matrix of 2-D embedding coordinates.}
#'   \item{\code{cell_metadata}}{Data frame of per-cell metadata extracted from
#'     the \code{SeuratObject}.}
#'   \item{\code{params}}{Named list of parameters used: \code{n_hvg},
#'     \code{knn_k}, \code{arrow_scale}, \code{point_colour}, \code{title}.}
#' }
#'
#' @name cpf
#' @aliases cpf-class
NULL


#' @export
print.cpf <- function(x, ...) {
  cat("Chromatin Potential Field <cpf>\n")
  cat("  Cells          :", nrow(x$embedding), "\n")
  cat("  HVGs used      :", length(x$hvg), "\n")
  cat("  kNN (CPF)      :", x$params$knn_k, "\n")
  cat("  Arrow scale    :", x$params$arrow_scale, "\n")
  cat("  Embedding dims :", paste(colnames(x$embedding), collapse = " / "), "\n")
  cat("  Colour by      :", x$params$point_colour, "\n")
  invisible(x)
}


#' @export
plot.cpf <- function(x, ...) {
  print(x$plot)
  invisible(x)
}


#' @export
summary.cpf <- function(object, ...) {
  mag <- object$cpf_df$magnitude
  cat("Chromatin Potential Field summary\n")
  cat("  Cells           :", nrow(object$embedding), "\n")
  cat("  HVGs            :", length(object$hvg), "\n")
  cat("  Arrow magnitude :",
      "min =", round(min(mag), 3),
      " median =", round(stats::median(mag), 3),
      " max =", round(max(mag), 3), "\n")
  invisible(object)
}


## ─────────────────────────────────────────────────────────────────────────────
## Main function
## ─────────────────────────────────────────────────────────────────────────────

#' Build the Chromatin Potential Field
#'
#' For every cell, identifies the \code{knn_k} cells whose RNA expression
#' profile is most correlated (Spearman) with that cell's scGAS vector.  The
#' vector from the cell to the centroid of those neighbours in 2-D embedding
#' space defines the Chromatin Potential Field (CPF) arrow.
#'
#' Returns a \code{\link{cpf}} object that holds the arrow data, the plot, the
#' HVG matrices, the embedding, and the parameters used — everything needed to
#' inspect or re-visualise the field without touching the original
#' \code{SeuratObject}.
#'
#' @param seurat_obj \code{SeuratObject} with \code{RNA} and \code{scGAS}
#'   assays (added by \code{\link{scgas_compute}}).
#' @param n_hvg Number of highly variable genes used to compute correlations.
#'   Default \code{2000}.
#' @param knn_k Number of RNA neighbours per cell.  Default \code{10}.
#' @param n_cores Cores for \code{parallel::mclapply}.  Default \code{1}.
#' @param arrow_scale Multiplicative scaling for arrow lengths.  Default
#'   \code{1}.
#' @param arrow_colour Colour for CPF arrows.  Default \code{"steelblue"}.
#' @param point_colour Metadata column name for background cell colouring.
#'   Default \code{"seurat_clusters"}.
#' @param point_size Size of background points.  Default \code{0.5}.
#' @param arrow_alpha Transparency of arrows.  Default \code{0.7}.
#' @param title Plot title.  Default \code{"Chromatin Potential Field"}.
#' @param embedding_coords Optional cells x 2 numeric matrix of 2-D
#'   coordinates.  \code{NULL} (default) uses \code{scgasumap} if present,
#'   else \code{umap}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return A \code{\link{cpf}} object with elements:
#'   \describe{
#'     \item{\code{cpf_df}}{Per-cell arrow data frame (x, y, dx, dy,
#'       magnitude, cell, optional metadata column).}
#'     \item{\code{plot}}{\code{ggplot2} object.}
#'     \item{\code{hvg}}{HVGs used.}
#'     \item{\code{scgas_hvg}}{genes x cells scGAS matrix for \code{hvg}.}
#'     \item{\code{rna_hvg}}{genes x cells RNA matrix for \code{hvg}.}
#'     \item{\code{embedding}}{cells x 2 embedding matrix.}
#'     \item{\code{cell_metadata}}{Data frame of cell metadata.}
#'     \item{\code{params}}{List of parameters used.}
#'   }
#'
#' @examples
#' \dontrun{
#' cpf_result <- scgas_chromatin_potential(
#'   seurat_obj   = obj,
#'   n_hvg        = 2000,
#'   knn_k        = 10,
#'   n_cores      = 8,
#'   arrow_scale  = 1.5,
#'   point_colour = "celltype"
#' )
#'
#' ## Print summary
#' print(cpf_result)
#'
#' ## Show the plot
#' plot(cpf_result)
#'
#' ## Subset to the strongest arrows
#' library(ggplot2)
#' strong <- cpf_result$cpf_df[
#'   cpf_result$cpf_df$magnitude > quantile(cpf_result$cpf_df$magnitude, 0.75), ]
#' }
#'
#' @importFrom Matrix t
#' @importFrom stats cor setNames var median
#' @importFrom parallel mclapply
#' @importFrom ggplot2 ggplot aes geom_point geom_segment theme_classic labs
#'   theme element_text
#' @importFrom grid arrow unit
#' @export
scgas_chromatin_potential <- function(seurat_obj,
                                     n_hvg            = 2000L,
                                     knn_k            = 10L,
                                     n_cores          = 1L,
                                     arrow_scale      = 1,
                                     arrow_colour     = "steelblue",
                                     point_colour     = "seurat_clusters",
                                     point_size       = 0.5,
                                     arrow_alpha      = 0.7,
                                     title            = "Chromatin Potential Field",
                                     embedding_coords = NULL,
                                     verbose          = TRUE) {

  stopifnot(inherits(seurat_obj, "Seurat"))

  ## ── Extract scGAS matrix ─────────────────────────────────────────────────
  if (!"scGAS" %in% names(seurat_obj@assays))
    stop("[scGAS] No 'scGAS' assay found. Run scgas_compute() first.")
  scgas_mat <- as.matrix(seurat_obj@assays$scGAS@data)

  ## ── Extract RNA matrix ───────────────────────────────────────────────────
  if (!"RNA" %in% names(seurat_obj@assays))
    stop("[scGAS] No RNA assay found in seurat_obj.")
  rna_mat <- as.matrix(seurat_obj@assays$RNA@data)

  ## ── Extract embedding ────────────────────────────────────────────────────
  if (is.null(embedding_coords)) {
    red_names <- names(seurat_obj@reductions)
    if ("scgasumap" %in% red_names) {
      embedding_coords <- seurat_obj@reductions$scgasumap@cell.embeddings[, 1:2]
    } else if ("umap" %in% red_names) {
      embedding_coords <- seurat_obj@reductions$umap@cell.embeddings[, 1:2]
    } else {
      stop("[scGAS] No UMAP reduction found. Supply embedding_coords.")
    }
  }
  colnames(embedding_coords) <- c("UMAP_1", "UMAP_2")

  ## ── Select HVGs ──────────────────────────────────────────────────────────
  if (verbose) message("[scGAS] Selecting HVGs for CPF computation")
  common_genes <- intersect(rownames(scgas_mat), rownames(rna_mat))
  if (length(common_genes) == 0)
    stop("[scGAS] No genes shared between scgas_mat and rna_mat.")

  scgas_vars <- apply(scgas_mat[common_genes, , drop = FALSE], 1, stats::var)
  rna_vars   <- apply(rna_mat[common_genes, , drop = FALSE],   1, stats::var)
  n_take     <- min(n_hvg, length(common_genes))
  hvg_scgas  <- names(sort(scgas_vars, decreasing = TRUE))[seq_len(n_take)]
  hvg_rna    <- names(sort(rna_vars,   decreasing = TRUE))[seq_len(n_take)]
  hvg        <- intersect(hvg_scgas, hvg_rna)

  if (length(hvg) < 10)
    warning("[scGAS] Only ", length(hvg), " HVGs shared. Results may be noisy.")
  if (verbose) message("[scGAS] Using ", length(hvg), " HVGs")

  gas_hvg <- t(scgas_mat[hvg, , drop = FALSE])   # cells x hvg
  rna_hvg <- t(rna_mat[hvg, , drop = FALSE])      # cells x hvg
  all_cells <- colnames(seurat_obj)

  ## ── Compute CPF arrows ───────────────────────────────────────────────────
  if (verbose) message("[scGAS] Computing CPF arrows (k = ", knn_k, ")")
  cpf_list <- parallel::mclapply(
    X        = as.list(all_cells),
    FUN      = .cpf_one_cell,
    gas_hvg  = gas_hvg,
    rna_hvg  = rna_hvg,
    embedding = embedding_coords,
    knn_k    = knn_k,
    mc.cores = n_cores
  )
  names(cpf_list) <- all_cells

  cpf_df <- do.call(rbind, lapply(all_cells, function(cell) {
    r <- cpf_list[[cell]]
    if (is.null(r)) return(NULL)
    data.frame(
      cell      = cell,
      x         = embedding_coords[cell, 1],
      y         = embedding_coords[cell, 2],
      dx        = r$dx * arrow_scale,
      dy        = r$dy * arrow_scale,
      magnitude = r$magnitude,
      stringsAsFactors = FALSE
    )
  }))
  rownames(cpf_df) <- cpf_df$cell

  if (point_colour %in% colnames(seurat_obj@meta.data))
    cpf_df[[point_colour]] <- as.character(
      seurat_obj@meta.data[cpf_df$cell, point_colour])

  ## ── Build plot ───────────────────────────────────────────────────────────
  bg_df <- data.frame(
    x = embedding_coords[all_cells, 1],
    y = embedding_coords[all_cells, 2]
  )
  if (point_colour %in% colnames(seurat_obj@meta.data))
    bg_df[[point_colour]] <- as.character(
      seurat_obj@meta.data[all_cells, point_colour])

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data    = bg_df,
      mapping = ggplot2::aes(x = x, y = y),
      colour  = "grey80",
      size    = point_size,
      alpha   = 0.6
    ) +
    ggplot2::geom_segment(
      data    = cpf_df,
      mapping = ggplot2::aes(x = x, y = y, xend = x + dx, yend = y + dy),
      arrow     = grid::arrow(length = grid::unit(0.08, "cm"), type = "closed"),
      alpha     = arrow_alpha,
      colour    = arrow_colour,
      linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text   = ggplot2::element_text(face = "bold"),
      axis.title  = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )

  ## ── Collect cell metadata ────────────────────────────────────────────────
  cell_metadata <- seurat_obj@meta.data[all_cells, , drop = FALSE]

  ## ── Construct and return cpf object ─────────────────────────────────────
  structure(
    list(
      cpf_df        = cpf_df,
      plot          = p,
      hvg           = hvg,
      scgas_hvg     = t(gas_hvg),   # back to genes x cells
      rna_hvg       = t(rna_hvg),   # back to genes x cells
      embedding     = embedding_coords,
      cell_metadata = cell_metadata,
      params        = list(
        n_hvg        = n_hvg,
        knn_k        = knn_k,
        arrow_scale  = arrow_scale,
        point_colour = point_colour,
        title        = title
      )
    ),
    class = "cpf"
  )
}


.cpf_one_cell <- function(cell, gas_hvg, rna_hvg, embedding, knn_k) {
  tryCatch({
    gas_vec    <- gas_hvg[cell, ]
    cors       <- apply(rna_hvg, 1, function(r) stats::cor(gas_vec, r, method = "spearman"))
    cors[cell] <- -Inf
    top_k      <- names(sort(cors, decreasing = TRUE))[seq_len(knn_k)]
    target_xy  <- colMeans(embedding[top_k, , drop = FALSE])
    self_xy    <- embedding[cell, ]
    delta      <- target_xy - self_xy
    list(dx = delta[1], dy = delta[2], magnitude = sqrt(sum(delta^2)))
  }, error = function(e) NULL)
}
