#' scGAS: Single-Cell Gene Activation Potential Inference
#'
#' @description
#' \pkg{scGAS} infers Gene Activation Potential (GAS) from single-cell
#' ATAC-seq data using a pre-built bulk-tissue reference map of cCRE-to-gene
#' associations and per-gene Lasso regression models.
#'
#' ## Main workflow
#'
#' Every function from Step 2 onwards takes a \code{SeuratObject} as its first
#' argument and returns the same object with results attached — either as cell
#' metadata, a new assay, or an entry in \code{@@misc}.  Step 5 is the only
#' exception: it returns a standalone \code{\link{cpf}} object.
#'
#' | Step | Function | Input → Output |
#' |------|----------|----------------|
#' | 1 | \code{\link{scgas_preprocess}} | Fragment files → \code{SeuratObject} with ATAC assay, LSI, UMAP, KNN graph |
#' | 2 | \code{\link{scgas_metacell}} | \code{SeuratObject} → adds \code{mc_membership} to cell metadata |
#' | 3 | \code{\link{scgas_train_models}} | \code{SeuratObject} → adds fitted models to \code{@@misc$scgas_models} |
#' | 4 | \code{\link{scgas_compute}} | \code{SeuratObject} → adds \code{scGAS} assay, \code{scgaspca}, \code{scgasumap} |
#' | 5 | \code{\link{scgas_chromatin_potential}} | \code{SeuratObject} → \code{\link{cpf}} object (arrows, plot, HVG matrices, embedding) |
#'
#' \code{\link{scgas_download_reference}} downloads the ENCODE reference files
#' from Zenodo into a local \code{data_dir} — run once after installation.
#' \code{\link{scgas_add_assay}} re-attaches a pre-computed scGAS matrix or
#' re-runs dimensionality reduction without recomputing GAS.
#'
#' @docType package
#' @name scGAS-package
#' @aliases scGAS
"_PACKAGE"
