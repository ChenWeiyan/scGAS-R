#' Train Per-Gene Lasso Models and Compute Metacell-Level GAS
#'
#' For each gene with enough cCRE-Gene associations in the dataset, trains a
#' Lasso regression model using ENCODE bulk DNase-seq and RNA-seq data, then
#' predicts Gene Activation Scores (GAS) at Metacell resolution.
#'
#' The fitted models and their associated metadata are stored in
#' \code{seurat_obj@@misc$scgas_models} and the updated object is returned.
#'
#' ## Automatic reference loading
#'
#' When \code{sig_associations}, \code{encode_rna_lognorm}, or
#' \code{encode_dnase_lognorm} are \code{NULL} (the default), the function
#' reads \code{seurat_obj@@misc$genome_version} and
#' \code{seurat_obj@@misc$data_dir} (set by \code{\link{scgas_preprocess}})
#' and loads the matching ENCODE reference files automatically.
#'
#' Supply explicit objects to override auto-loading for any of the three.
#'
#' @param seurat_obj A preprocessed \code{SeuratObject} from
#'   \code{\link{scgas_preprocess}} with Metacell membership in
#'   \code{seurat_obj$mc_membership} (added by \code{\link{scgas_metacell}}).
#' @param sig_associations Named list of per-gene association objects.  Each
#'   element must contain a \code{Sig.df} data frame with a \code{CRE.sig}
#'   column.  \code{NULL} (default) auto-loads from \code{data_dir}.
#' @param encode_rna_lognorm Matrix (genes x samples) of log-normalised
#'   bulk RNA-seq values from ENCODE tissues.  \code{NULL} auto-loads.
#' @param encode_dnase_lognorm Matrix (cCREs x samples) of log-normalised
#'   bulk DNase-seq values.  \code{NULL} auto-loads.
#' @param cre_bed_path Path to the reference cCRE BED file.  \code{NULL}
#'   uses the \code{GRanges} cached in \code{seurat_obj@@misc$ref_cre_gr}.
#' @param encode_lib_size Total cCRE count in the ENCODE reference used to
#'   compute the library-size scaling factor.  Default \code{2081250L}.
#' @param min_associations Minimum cCRE-Gene associations required to model a
#'   gene.  \code{NULL} (default) uses the median across all genes.
#' @param lasso_alpha Elastic-net mixing parameter.  \code{1} = pure Lasso,
#'   \code{0} = Ridge.  Default \code{1}.
#' @param deviance_quantile Lambda selection threshold.  Default \code{0.95}.
#' @param n_cores Cores for \code{parallel::mclapply}.  Default \code{1}.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return The input \code{SeuratObject} with fitted models stored in
#'   \code{seurat_obj@@misc$scgas_models}.  Each element is named by gene and
#'   contains \code{glmnetfit}, \code{lambda_choose}, \code{GAS}, and
#'   \code{cCRE_used}.
#'
#' @examples
#' \dontrun{
#' ## Simplest call — reference files loaded automatically
#' obj <- scgas_train_models(obj, n_cores = 8)
#'
#' ## Override one reference explicitly
#' obj <- scgas_train_models(obj, sig_associations = my_sig_ass, n_cores = 8)
#'
#' message(length(obj@misc$scgas_models), " genes modelled")
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges findOverlaps sort seqnames start end
#' @importFrom S4Vectors queryHits
#' @importFrom Matrix rowSums
#' @importFrom MatrixGenerics rowSums2
#' @importFrom parallel mclapply
#' @importFrom glmnetUtils glmnet
#' @importFrom stats predict median
#' @export
scgas_train_models <- function(seurat_obj,
                               sig_associations     = NULL,
                               encode_rna_lognorm   = NULL,
                               encode_dnase_lognorm = NULL,
                               cre_bed_path         = NULL,
                               encode_lib_size      = 2081250L,
                               min_associations     = NULL,
                               lasso_alpha          = 1,
                               deviance_quantile    = 0.95,
                               n_cores              = 1L,
                               verbose              = TRUE) {

  stopifnot(
    inherits(seurat_obj, "Seurat"),
    is.numeric(encode_lib_size), encode_lib_size > 0,
    is.numeric(lasso_alpha), lasso_alpha >= 0, lasso_alpha <= 1,
    is.numeric(deviance_quantile), deviance_quantile > 0, deviance_quantile <= 1
  )

  ## ── Auto-load reference files if not supplied ─────────────────────────────
  genome_version <- seurat_obj@misc$genome_version
  data_dir       <- seurat_obj@misc$data_dir

  if (is.null(sig_associations)) {
    sig_associations <- .autoload_reference(
      "sig_ass", genome_version, data_dir, verbose)
  }
  if (is.null(encode_rna_lognorm)) {
    encode_rna_lognorm <- .autoload_reference(
      "rna_enc", genome_version, data_dir, verbose)
  }
  if (is.null(encode_dnase_lognorm)) {
    encode_dnase_lognorm <- .autoload_reference(
      "dnase_enc", genome_version, data_dir, verbose)
  }

  stopifnot(
    is.list(sig_associations), !is.null(names(sig_associations)),
    is.matrix(encode_rna_lognorm) || is.data.frame(encode_rna_lognorm),
    is.matrix(encode_dnase_lognorm)
  )

  ## ── Resolve Metacell labels ───────────────────────────────────────────────
  if ("mc_membership" %in% colnames(seurat_obj@meta.data)) {
    cl_labels <- seurat_obj$mc_membership
  } else {
    stop("[scGAS] 'mc_membership' not found in metadata. ",
         "Run scgas_metacell() first.")
  }

  if (verbose) message("[scGAS] Matching cCRE names to dataset features")
  ref_cre_gr <- .load_ref_cre(seurat_obj, cre_bed_path, verbose)

  obj_ranges <- seurat_obj@assays$ATAC@ranges
  ov <- GenomicRanges::findOverlaps(ref_cre_gr, obj_ranges)
  obj_ranges$name <- ref_cre_gr[S4Vectors::queryHits(ov)]$name
  names(obj_ranges) <- obj_ranges$name

  all_ass_cre  <- unique(unlist(lapply(sig_associations, function(x) x$Sig.df$CRE.sig)))
  all_ass_gr   <- GenomicRanges::sort(obj_ranges[which(obj_ranges$name %in% all_ass_cre)])
  all_ass_names <- paste(
    as.character(GenomicRanges::seqnames(all_ass_gr)),
    GenomicRanges::start(all_ass_gr),
    GenomicRanges::end(all_ass_gr),
    sep = "-"
  )

  if (verbose) message("[scGAS] Aggregating counts into Metacells")
  cl_levels <- levels(as.factor(cl_labels))

  count_mat  <- seurat_obj@assays$ATAC@counts
  feat_inter <- intersect(rownames(count_mat), all_ass_names)
  count_mat  <- count_mat[feat_inter, , drop = FALSE]

  atac_cl_mat <- matrix(
    0L,
    nrow = length(feat_inter),
    ncol = length(cl_levels),
    dimnames = list(feat_inter, paste0("Cluster-", cl_levels))
  )
  for (cl in cl_levels) {
    sel <- which(cl_labels == cl)
    atac_cl_mat[, paste0("Cluster-", cl)] <-
      MatrixGenerics::rowSums2(count_mat[, sel, drop = FALSE])
  }

  sf <- floor(encode_lib_size * length(feat_inter) / nrow(encode_dnase_lognorm))
  if (verbose) message("[scGAS] Library-size scaling factor: ", sf)
  atac_cl_lognorm <- log2(t(t(atac_cl_mat) / colSums(atac_cl_mat)) * sf + 1)

  name_map <- stats::setNames(all_ass_gr$name, all_ass_names)
  rownames(atac_cl_lognorm) <- name_map[rownames(atac_cl_lognorm)]

  ass_counts <- vapply(sig_associations, function(x) {
    length(intersect(all_ass_gr$name, x$Sig.df$CRE.sig))
  }, integer(1L))
  thr       <- if (is.null(min_associations)) stats::median(ass_counts) else as.integer(min_associations)
  sel_genes <- names(which(ass_counts > thr))

  if (length(sel_genes) == 0)
    stop("[scGAS] No genes passed the association threshold (", thr,
         "). Try lowering min_associations.")
  if (verbose) message("[scGAS] Training models for ", length(sel_genes),
                       " genes (threshold = ", thr, " associations)")

  trained_models <- parallel::mclapply(
    X                    = as.list(sel_genes),
    FUN                  = .train_one_gene,
    sig_associations     = sig_associations,
    all_ass_gr           = all_ass_gr,
    encode_rna_lognorm   = encode_rna_lognorm,
    encode_dnase_lognorm = encode_dnase_lognorm,
    atac_cl_lognorm      = atac_cl_lognorm,
    lasso_alpha          = lasso_alpha,
    deviance_quantile    = deviance_quantile,
    mc.cores             = n_cores
  )
  names(trained_models) <- sel_genes

  failed <- vapply(trained_models, is.null, logical(1L))
  if (any(failed)) {
    warning("[scGAS] Training failed for ", sum(failed), " gene(s): ",
            paste(sel_genes[failed], collapse = ", "))
    trained_models <- trained_models[!failed]
  }

  attr(trained_models, "ATAC_CL_lognorm") <- atac_cl_lognorm
  attr(trained_models, "cl_levels")       <- cl_levels
  attr(trained_models, "all_ass_gr")      <- all_ass_gr

  if (verbose) message("[scGAS] Model training complete: ",
                       length(trained_models), " genes.")

  ## ── Store in SeuratObject and return ────────────────────────────────────
  seurat_obj@misc$scgas_models <- trained_models
  return(seurat_obj)
}


# --------------------------------------------------------------------------
# .autoload_reference: load an ENCODE reference file from data_dir.
# Errors with a clear message if the file cannot be found.
# --------------------------------------------------------------------------
#' @keywords internal
.autoload_reference <- function(type, genome_version, data_dir, verbose) {
  if (is.null(data_dir))
    stop("[scGAS] Cannot auto-load '", type, "': data_dir not set. ",
         "Run scgas_preprocess() first, or pass the reference object explicitly.")
  if (is.null(genome_version))
    stop("[scGAS] Cannot auto-load '", type, "': genome_version not set. ",
         "Run scgas_preprocess() first, or pass the reference object explicitly.")

  path <- .find_reference_file(type, genome_version, data_dir)
  if (is.null(path))
    stop("[scGAS] Cannot find reference file '", type, "' for genome '",
         genome_version, "' in data_dir '", data_dir, "'.\n",
         "Download the matching ENCODE reference files and place them in data_dir.")

  if (verbose) message("[scGAS] Loading ", type, ": ", basename(path))
  readRDS(path)
}


.load_ref_cre <- function(seurat_obj, cre_bed_path, verbose) {
  if (!is.null(seurat_obj@misc$ref_cre_gr)) {
    if (verbose) message("[scGAS] Using cached ref_cre_gr from seurat_obj@misc")
    return(seurat_obj@misc$ref_cre_gr)
  }
  if (!is.null(cre_bed_path)) {
    if (!file.exists(cre_bed_path))
      stop("[scGAS] cre_bed_path not found: ", cre_bed_path)
    if (verbose) message("[scGAS] Loading reference cCRE BED: ", cre_bed_path)
    gr <- rtracklayer::import(cre_bed_path)
    names(gr) <- gr$name
    return(gr)
  }
  stop("[scGAS] Provide cre_bed_path or run scgas_preprocess() first ",
       "(which caches ref_cre_gr in seurat_obj@misc).")
}


.train_one_gene <- function(gene_name,
                            sig_associations,
                            all_ass_gr,
                            encode_rna_lognorm,
                            encode_dnase_lognorm,
                            atac_cl_lognorm,
                            lasso_alpha,
                            deviance_quantile) {
  tryCatch({
    gene_cres <- sig_associations[[gene_name]]$Sig.df$CRE.sig
    cre_v     <- intersect(all_ass_gr$name, gene_cres)
    if (length(cre_v) == 0) return(NULL)

    train_df <- data.frame(
      RNA = encode_rna_lognorm[gene_name, ],
      t(encode_dnase_lognorm[cre_v, , drop = FALSE])
    )

    fit <- glmnetUtils::glmnet(
      RNA ~ .,
      data   = train_df,
      family = "gaussian",
      alpha  = lasso_alpha
    )
    class(fit) <- "glmnet"

    which_max     <- which(fit$dev.ratio >= max(fit$dev.ratio) * deviance_quantile)[1]
    lambda_choose <- fit$lambda[which_max]

    newx <- t(atac_cl_lognorm[cre_v, , drop = FALSE])
    gas  <- as.vector(stats::predict(fit, newx = newx, type = "response",
                                     s = lambda_choose))
    gas  <- gas - min(gas)

    list(glmnetfit = fit, lambda_choose = lambda_choose,
         GAS = gas, cCRE_used = cre_v)
  }, error = function(e) NULL)
}
