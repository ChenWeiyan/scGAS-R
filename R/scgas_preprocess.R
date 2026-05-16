#' Preprocess scATAC-seq Fragment Files into a SeuratObject
#'
#' Builds a \code{SeuratObject} with an \code{ATAC} assay whose features are
#' the reference cCREs, computes QC metrics, runs TF-IDF/LSI dimensionality
#' reduction, and builds a KNN graph ready for \code{\link{scgas_metacell}}.
#'
#' ## Automatic genome and reference resolution
#'
#' When \code{genome_version = NULL} (the default), the function reads the
#' \code{#} comment header of the fragment file (written by CellRanger /
#' CellRanger-ARC) and extracts the genome assembly token (e.g.
#' \code{GRCh38}, \code{GRCh37}).  It then maps that token to the
#' corresponding UCSC name (\code{"hg38"}, \code{"hg19"}, \ldots) and looks
#' for matching reference files in \code{data_dir}.
#'
#' If the genome cannot be detected, a warning is issued and \code{"hg19"} is
#' used as the fallback.
#'
#' When \code{genome_version} is set explicitly, detection is skipped
#' entirely.
#'
#' The detected or supplied genome and the resolved \code{data_dir} are stored
#' in \code{seurat_obj@@misc} so that \code{\link{scgas_train_models}} can
#' load the correct ENCODE reference files automatically.
#'
#' ## Annotation
#'
#' Gene annotations are resolved automatically from the detected genome
#' assembly using the appropriate \code{EnsDb} package
#' (\code{EnsDb.Hsapiens.v75} for hg19, \code{EnsDb.Hsapiens.v86} for hg38).
#' \code{seqlevelsStyle} and \code{genome} are set to match the fragment file.
#' No manual annotation step is needed.
#'
#' To override, pass a \code{GRanges} object via the \code{annotation}
#' argument; style and genome will still be adjusted automatically.
#'
#' @param fragment_path Path to a fragment file (\code{.tsv.gz}).  For
#'   multi-replicate datasets supply a character vector — one element per
#'   replicate.  A single path is the most common case.
#' @param cre_bed_path Path to the reference cCRE BED file.  \code{NULL}
#'   (default) looks for \code{{genome}_500bp_CRE.bed} inside \code{data_dir}.
#' @param data_dir Path to the directory containing scGAS reference files
#'   (\code{hg19_500bp_CRE.bed}, ENCODE \code{.gz} files, etc.).  \code{NULL}
#'   (default) searches for \code{"data"} and \code{"../data"} relative to
#'   the current working directory.
#' @param genome_version Genome assembly string (\code{"hg19"}, \code{"hg38"},
#'   \code{"mm10"}, \ldots).  \code{NULL} (default) auto-detects from the
#'   fragment file header.  Set explicitly to skip detection.
#' @param rna_matrix Optional sparse RNA count matrix (genes x cells) for
#'   multiome data.
#' @param meta_data Optional data frame of per-cell metadata (row names =
#'   cell barcodes).
#' @param annotation Optional \code{GRanges} of gene annotations.  When
#'   \code{NULL} (default) the appropriate \code{EnsDb} package is loaded
#'   automatically based on the detected genome.  \code{seqlevelsStyle} and
#'   \code{genome} are always set automatically.
#' @param cell_barcodes Optional character vector of barcodes to retain.
#' @param replicate_ids Character vector of replicate labels when
#'   \code{fragment_path} is a vector.  Defaults to \code{"rep1"},
#'   \code{"rep2"}, \ldots
#' @param min_cells_per_feature Minimum cells a cCRE must be detected in.
#'   Default \code{1}.
#' @param min_features_per_cell Minimum cCREs detected per cell.  Default
#'   \code{1}.
#' @param lsi_dims LSI components for UMAP and neighbour graph.  Default
#'   \code{2:30}.
#' @param umap_n_components Number of UMAP components.  Default \code{3}.
#' @param n_cores Cores for \code{parallel::mclapply}.  Default \code{1}.
#' @param verbose Logical.  Default \code{TRUE}.
#' @param out_dir Path to a base output directory shared by all pipeline steps.
#'   Stored in \code{seurat_obj@@misc$out_dir} so downstream functions
#'   (\code{scgas_metacell}, \code{scgas_train_models}, \code{scgas_compute})
#'   inherit it automatically.  \code{NULL} disables disk caching for all
#'   steps.
#' @param run_name Character label for this analysis run.  All intermediate
#'   files are written under \code{out_dir/run_name/}.  Stored in
#'   \code{seurat_obj@@misc$run_name}.  Default \code{"scgas_run"}.
#' @param save_obj Logical.  Save the returned \code{SeuratObject} to
#'   \code{out_dir/run_name/seurat_preprocessed.rds}.  Requires \code{out_dir}
#'   to be set.  Default \code{FALSE}.
#'
#' @return A \code{SeuratObject} with:
#'   \describe{
#'     \item{\code{ATAC}}{Chromatin assay with raw cCRE counts and LSI.}
#'     \item{\code{RNA}}{(Optional) RNA assay with log-normalised data and PCA.}
#'     \item{\code{umap}}{UMAP embedding.}
#'     \item{\code{ATAC_snn}}{KNN graph for \code{\link{scgas_metacell}}.}
#'     \item{\code{misc$ref_cre_gr}}{Reference cCRE \code{GRanges}.}
#'     \item{\code{misc$genome_version}}{Resolved genome assembly string.}
#'     \item{\code{misc$data_dir}}{Resolved path to reference data directory.}
#'     \item{\code{misc$out_dir}}{Base output directory (or \code{NULL}).}
#'     \item{\code{misc$run_name}}{Run label used by downstream functions.}
#'   }
#'
#' @examples
#' \dontrun{
#' ## Genome is auto-detected; annotation is loaded automatically from EnsDb
#' obj <- scgas_preprocess(
#'   fragment_path = "datasets/10k_Human_Brain_MO_gemx_atac_fragments.tsv.gz",
#'   data_dir      = "data"
#' )
#'
#' ## Override genome detection explicitly (annotation still auto-loaded)
#' obj <- scgas_preprocess(
#'   fragment_path  = "datasets/my_atac_fragments.tsv.gz",
#'   data_dir       = "data",
#'   genome_version = "hg38"
#' )
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom Signac CreateFragmentObject FeatureMatrix CreateChromatinAssay
#'   NucleosomeSignal TSSEnrichment FindTopFeatures RunTFIDF RunSVD
#' @importFrom Seurat CreateSeuratObject CreateAssayObject FindNeighbors
#'   RunUMAP AddMetaData
#' @importFrom Matrix Matrix rowSums
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomeInfoDb seqlevels
#' @export
scgas_preprocess <- function(fragment_path,
                             cre_bed_path           = NULL,
                             data_dir               = NULL,
                             genome_version         = NULL,
                             rna_matrix             = NULL,
                             meta_data              = NULL,
                             annotation             = NULL,
                             cell_barcodes          = NULL,
                             replicate_ids          = NULL,
                             min_cells_per_feature  = 1L,
                             min_features_per_cell  = 1L,
                             lsi_dims               = 2:30,
                             umap_n_components      = 3L,
                             n_cores                = 1L,
                             verbose                = TRUE,
                             out_dir                = NULL,
                             run_name               = "scgas_run",
                             save_obj               = FALSE) {

  ## ── Validate fragment paths ───────────────────────────────────────────────
  stopifnot(is.character(fragment_path), length(fragment_path) >= 1)
  for (fp in fragment_path) {
    if (!file.exists(fp)) stop("[scGAS] Fragment file not found: ", fp)
  }

  n_reps <- length(fragment_path)
  if (is.null(replicate_ids))
    replicate_ids <- if (n_reps == 1) "rep1" else paste0("rep", seq_len(n_reps))
  stopifnot(length(replicate_ids) == n_reps)

  ## ── Resolve genome version ────────────────────────────────────────────────
  if (!is.null(genome_version)) {
    if (verbose) message("[scGAS] Using user-specified genome: ", genome_version)
  } else {
    genome_version <- .detect_genome_from_fragments(fragment_path[1], verbose)
    if (is.null(genome_version)) {
      warning("[scGAS] Could not detect genome assembly from fragment file header. ",
              "Defaulting to 'hg19'. Set genome_version explicitly to suppress this warning.")
      genome_version <- "hg19"
    }
  }

  ## ── Auto-load annotation if not supplied ─────────────────────────────────
  if (is.null(annotation))
    annotation <- .auto_annotation(genome_version, verbose)

  ## ── Resolve data directory ────────────────────────────────────────────────
  data_dir <- .resolve_data_dir(data_dir, verbose)

  ## ── Resolve CRE BED path ─────────────────────────────────────────────────
  if (is.null(cre_bed_path)) {
    cre_bed_path <- .find_reference_file("cre_bed", genome_version, data_dir)
    if (is.null(cre_bed_path))
      stop("[scGAS] Cannot find CRE BED file for genome '", genome_version,
           "' in data_dir '", data_dir, "'.\n",
           "Expected: ", file.path(data_dir, .cre_bed_name(genome_version)))
    if (verbose) message("[scGAS] Using CRE BED: ", cre_bed_path)
  }
  stopifnot(file.exists(cre_bed_path))

  ## ── Load reference BED ────────────────────────────────────────────────────
  if (verbose) message("[scGAS] Loading reference cCRE BED")
  ref_cre_gr <- rtracklayer::import(cre_bed_path)
  names(ref_cre_gr) <- ref_cre_gr$name

  ## ── Detect fragment chromosome style and harmonise BED ───────────────────
  frag_style <- .detect_chr_style(fragment_path[1], verbose)
  ref_cre_gr <- .harmonise_seqlevels(ref_cre_gr, frag_style, verbose)

  ## ── Prepare annotation ────────────────────────────────────────────────────
  if (!is.null(annotation))
    annotation <- .prepare_annotation(annotation, frag_style, genome_version, verbose)

  ## ── Build per-replicate count matrices ────────────────────────────────────
  count_list <- vector("list", n_reps)
  frag_list  <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    if (verbose) message("[scGAS] Processing replicate ", replicate_ids[i])
    frag_list[[i]] <- Signac::CreateFragmentObject(fragment_path[i])
    cm <- Signac::FeatureMatrix(fragments = frag_list[[i]], features = ref_cre_gr)
    if (!is.null(cell_barcodes))
      cm <- cm[, intersect(colnames(cm), cell_barcodes), drop = FALSE]
    count_list[[i]] <- cm
  }

  merged_counts <- if (n_reps == 1) count_list[[1]] else do.call(cbind, count_list)

  zero_rows <- which(Matrix::rowSums(merged_counts) == 0)
  if (length(zero_rows) > 0) {
    if (verbose) message("[scGAS] Removing ", length(zero_rows), " all-zero cCREs")
    merged_counts <- merged_counts[-zero_rows, , drop = FALSE]
    ref_cre_gr    <- ref_cre_gr[-zero_rows]
  }

  ## ── Build per-replicate SeuratObjects ─────────────────────────────────────
  seurat_list <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    rep_counts  <- count_list[[i]][rownames(merged_counts), , drop = FALSE]
    chrom_assay <- Signac::CreateChromatinAssay(
      fragments    = frag_list[[i]],
      counts       = rep_counts,
      ranges       = ref_cre_gr,
      min.cells    = min_cells_per_feature,
      min.features = min_features_per_cell,
      annotation   = annotation
    )
    seurat_list[[i]] <- Seurat::CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
  }

  obj <- if (n_reps == 1) {
    seurat_list[[1]]
  } else {
    merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = replicate_ids)
  }

  ## ── Optional RNA assay ────────────────────────────────────────────────────
  if (!is.null(rna_matrix)) {
    if (verbose) message("[scGAS] Adding RNA assay")
    rna_matrix <- Matrix::Matrix(rna_matrix, sparse = TRUE)
    rna_matrix <- rna_matrix[, intersect(colnames(rna_matrix), colnames(obj)), drop = FALSE]
    obj[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_matrix)
    Seurat::DefaultAssay(obj) <- "RNA"
    obj <- Seurat::FindVariableFeatures(obj, verbose = verbose)
    obj <- Seurat::ScaleData(obj, verbose = verbose)
    obj <- Seurat::RunPCA(obj, verbose = verbose)
  }

  ## ── Optional per-cell metadata ────────────────────────────────────────────
  if (!is.null(meta_data)) {
    shared <- intersect(rownames(meta_data), colnames(obj))
    obj <- Seurat::AddMetaData(obj, metadata = meta_data[shared, , drop = FALSE])
  }

  ## ── QC, TF-IDF/LSI, UMAP, KNN ────────────────────────────────────────────
  Seurat::DefaultAssay(obj) <- "ATAC"
  if (verbose) message("[scGAS] Computing QC metrics")
  obj <- Signac::NucleosomeSignal(obj)
  obj <- Signac::TSSEnrichment(obj)

  if (verbose) message("[scGAS] Running LSI")
  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q10", verbose = verbose)
  obj <- Signac::RunTFIDF(obj, verbose = verbose)
  obj <- Signac::RunSVD(obj, verbose = verbose)

  if (verbose) message("[scGAS] Running UMAP (", umap_n_components, " components)")
  obj <- Seurat::RunUMAP(
    object       = obj,
    reduction    = "lsi",
    dims         = lsi_dims,
    n.components = as.integer(umap_n_components),
    verbose      = verbose
  )

  if (verbose) message("[scGAS] Building KNN graph")
  obj <- Seurat::FindNeighbors(obj, reduction = "lsi", dims = lsi_dims,
                               verbose = FALSE)

  ## ── Cache reference info in misc ─────────────────────────────────────────
  obj@misc$ref_cre_gr      <- ref_cre_gr
  obj@misc$genome_version  <- genome_version
  obj@misc$data_dir        <- data_dir
  obj@misc$out_dir         <- out_dir
  obj@misc$run_name        <- run_name

  if (verbose) message("[scGAS] Preprocessing complete: ",
                       ncol(obj), " cells, ", nrow(obj), " cCREs  [", genome_version, "]")

  if (save_obj && !is.null(out_dir)) {
    run_dir  <- file.path(out_dir, run_name)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    obj_path <- file.path(run_dir, "seurat_preprocessed.rds")
    if (verbose) message("[scGAS] Saving preprocessed object to ", obj_path)
    saveRDS(obj, obj_path)
  }

  return(obj)
}


## ─────────────────────────────────────────────────────────────────────────────
## Reference file naming conventions
## ─────────────────────────────────────────────────────────────────────────────

## CRE BED: {genome}_500bp_CRE.bed
.cre_bed_name <- function(genome) sprintf("%s_500bp_CRE.bed", genome)

## Significance association list: genome-specific first, then generic fallback
.sigass_names <- function(genome) {
  c(sprintf("ENCODE_167Tissue_SigAss_%s_pVal.05.gz",  genome),
    sprintf("ENCODE_167Tissue_SigAss_%s_pVal.05.rds", genome),
    "ENCODE_167Tissue_SigAss_pVal.05.gz",
    "ENCODE_167Tissue_SigAss_pVal.05.rds")
}

## ENCODE RNA: genome-agnostic (same file for all assemblies)
.rna_enc_names <- function(genome) {   # nolint: unused arg kept for API consistency
  c("ENCODE_167Tissue_RNA_LogNorm.gz",
    "ENCODE_167Tissue_RNA_LogNorm.rds")
}

## ENCODE DNase: genome-specific first, then generic fallback
.dnase_enc_names <- function(genome) {
  c(sprintf("ENCODE_167Tissue_CRE_%s_LogNorm.gz",  genome),
    sprintf("ENCODE_167Tissue_CRE_%s_LogNorm.rds", genome),
    "ENCODE_167Tissue_CRE_LogNorm.gz",
    "ENCODE_167Tissue_CRE_LogNorm.rds")
}

## Unified lookup ---------------------------------------------------------
#' @keywords internal
.find_reference_file <- function(type, genome, data_dir) {
  candidates <- switch(type,
    cre_bed   = .cre_bed_name(genome),
    sig_ass   = .sigass_names(genome),
    rna_enc   = .rna_enc_names(genome),
    dnase_enc = .dnase_enc_names(genome),
    stop("[scGAS] Unknown reference type: ", type)
  )
  for (nm in candidates) {
    p <- file.path(data_dir, nm)
    if (file.exists(p)) return(p)
  }
  return(NULL)
}


## ─────────────────────────────────────────────────────────────────────────────
## Internal helpers
## ─────────────────────────────────────────────────────────────────────────────

# --------------------------------------------------------------------------
# .resolve_data_dir: find the reference data directory.
# Tries explicit path first, then searches common relative locations.
# --------------------------------------------------------------------------
#' @keywords internal
.resolve_data_dir <- function(data_dir, verbose = TRUE) {
  if (!is.null(data_dir)) {
    if (!dir.exists(data_dir))
      stop("[scGAS] data_dir not found: ", data_dir)
    return(normalizePath(data_dir))
  }
  for (candidate in c("data", "../data")) {
    if (dir.exists(candidate)) {
      found <- normalizePath(candidate)
      if (verbose) message("[scGAS] Found data directory: ", found)
      return(found)
    }
  }
  stop("[scGAS] Cannot locate reference data directory. ",
       "Set data_dir = '/path/to/scGAS/data'.")
}


# --------------------------------------------------------------------------
# .detect_genome_from_fragments: parse the comment header of a fragment file
# (written by CellRanger / CellRanger-ARC) for a genome assembly token.
# Returns a UCSC-style string ("hg38", "hg19", ...) or NULL.
# --------------------------------------------------------------------------
#' @keywords internal
.detect_genome_from_fragments <- function(fragment_path, verbose = TRUE) {
  ## Map of patterns (case-insensitive) found in fragment headers to UCSC names
  pattern_map <- c(
    "GRCh38" = "hg38",  "hg38"   = "hg38",
    "GRCh37" = "hg19",  "hg19"   = "hg19",
    "GRCm39" = "mm39",  "mm39"   = "mm39",
    "GRCm38" = "mm10",  "mm10"   = "mm10",
    "GRCm37" = "mm9",   "mm9"    = "mm9",
    "Rnor_7" = "rn7",   "rn7"    = "rn7",
    "Rnor_6" = "rn6",   "rn6"    = "rn6",
    "dm6"    = "dm6",
    "ce11"   = "ce11"
  )

  tryCatch({
    con <- if (grepl("\\.gz$", fragment_path, ignore.case = TRUE))
      gzcon(file(fragment_path, "rb"))
    else
      file(fragment_path, "r")
    on.exit(close(con), add = TRUE)

    header_lines <- character(0)
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0 || !startsWith(line, "#")) break
      header_lines <- c(header_lines, line)
    }

    if (length(header_lines) == 0) {
      if (verbose) message("[scGAS] No header lines found in fragment file.")
      return(NULL)
    }

    all_text <- paste(header_lines, collapse = " ")
    for (pattern in names(pattern_map)) {
      if (grepl(pattern, all_text, ignore.case = TRUE)) {
        genome <- pattern_map[[pattern]]
        if (verbose) message("[scGAS] Detected genome from fragment header: ",
                             genome, "  (matched '", pattern, "')")
        return(genome)
      }
    }

    if (verbose) message("[scGAS] Fragment header found but no known genome token detected.")
    return(NULL)

  }, error = function(e) {
    if (verbose) message("[scGAS] Could not read fragment header: ", conditionMessage(e))
    return(NULL)
  })
}


# --------------------------------------------------------------------------
# .detect_chr_style: read the first data line and return "UCSC", "NCBI",
# or "unknown".
# --------------------------------------------------------------------------
#' @keywords internal
.detect_chr_style <- function(fragment_path, verbose = TRUE) {
  style <- tryCatch({
    con <- if (grepl("\\.gz$", fragment_path, ignore.case = TRUE))
      gzcon(file(fragment_path, "rb"))
    else
      file(fragment_path, "r")
    on.exit(close(con), add = TRUE)
    first_chrom <- NULL
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0) break
      if (!startsWith(line, "#")) {
        first_chrom <- strsplit(line, "\t")[[1]][1]
        break
      }
    }
    if (is.null(first_chrom)) "unknown"
    else if (startsWith(first_chrom, "chr")) "UCSC"
    else "NCBI"
  }, error = function(e) "unknown")

  if (style == "unknown" && verbose)
    message("[scGAS] Could not detect chromosome style from fragment file.")
  else if (verbose)
    message("[scGAS] Fragment file chromosome style: ", style)

  return(style)
}


# --------------------------------------------------------------------------
# .harmonise_seqlevels: adjust BED GRanges chr naming to match fragment style.
# --------------------------------------------------------------------------
#' @keywords internal
.harmonise_seqlevels <- function(gr, frag_style, verbose = TRUE) {
  if (frag_style == "unknown") {
    if (verbose) message("[scGAS] Leaving BED seqlevels unchanged (style unknown).")
    return(gr)
  }

  bed_style <- if (any(startsWith(as.character(GenomicRanges::seqnames(gr)), "chr")))
    "UCSC" else "NCBI"

  if (bed_style == frag_style) {
    if (verbose) message("[scGAS] BED chromosome style consistent (", frag_style, ").")
    return(gr)
  }

  if (frag_style == "NCBI") {
    if (verbose) message("[scGAS] Stripping 'chr' prefix from BED seqlevels.")
    GenomeInfoDb::seqlevels(gr) <- gsub("^chr", "", GenomeInfoDb::seqlevels(gr))
  } else {
    if (verbose) message("[scGAS] Adding 'chr' prefix to BED seqlevels.")
    GenomeInfoDb::seqlevels(gr) <- paste0("chr", GenomeInfoDb::seqlevels(gr))
  }
  return(gr)
}


# --------------------------------------------------------------------------
# .auto_annotation: load the appropriate EnsDb package for a given genome
# and return a GRanges annotation object. Returns NULL (with a warning) if
# the required package is not installed.
# --------------------------------------------------------------------------
#' @keywords internal
.auto_annotation <- function(genome_version, verbose = TRUE) {
  pkg_map <- c(
    hg19 = "EnsDb.Hsapiens.v75",
    hg38 = "EnsDb.Hsapiens.v86",
    mm10 = "EnsDb.Mmusculus.v79",
    mm39 = "EnsDb.Mmusculus.v109"
  )
  pkg <- pkg_map[[genome_version]]
  if (is.null(pkg)) {
    if (verbose)
      message("[scGAS] No EnsDb package mapped for genome '", genome_version,
              "'. Proceeding without annotation (TSS enrichment unavailable).")
    return(NULL)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning("[scGAS] Package '", pkg, "' is required for automatic annotation ",
            "with genome '", genome_version, "' but is not installed.\n",
            "Install it with: BiocManager::install(\"", pkg, "\")\n",
            "Proceeding without annotation (TSS enrichment unavailable).")
    return(NULL)
  }
  if (verbose) message("[scGAS] Auto-loading annotation from ", pkg)
  db  <- get(pkg, envir = asNamespace(pkg))
  ann <- Signac::GetGRangesFromEnsDb(ensdb = db, verbose = FALSE)
  return(ann)
}


# --------------------------------------------------------------------------
# .resolve_run_dir: resolve out_dir / run_name from the function parameter
# (highest priority) or from seurat_obj@misc (inherited from a prior step).
# Returns a list with out_dir, run_name, run_dir, and use_cache flag.
# --------------------------------------------------------------------------
#' @keywords internal
.resolve_run_dir <- function(out_dir, run_name, misc, verbose) {
  resolved_out  <- if (!is.null(out_dir))  out_dir  else misc$out_dir
  resolved_name <- if (!is.null(run_name)) run_name else misc$run_name

  if (is.null(resolved_out))
    return(list(out_dir = NULL, run_name = NULL, run_dir = NULL, use_cache = FALSE))

  if (is.null(resolved_name)) resolved_name <- "scgas_run"

  if (verbose && is.null(out_dir) && !is.null(misc$out_dir))
    message("[scGAS] Using run: ", file.path(resolved_out, resolved_name))

  list(
    out_dir   = resolved_out,
    run_name  = resolved_name,
    run_dir   = file.path(resolved_out, resolved_name),
    use_cache = TRUE
  )
}


# --------------------------------------------------------------------------
# .prepare_annotation: auto-set seqlevels style and genome on annotation.
# --------------------------------------------------------------------------
#' @keywords internal
.prepare_annotation <- function(annotation, frag_style, genome_str, verbose = TRUE) {
  if (frag_style != "unknown") {
    ann_style <- if (any(startsWith(
      as.character(GenomicRanges::seqnames(annotation)), "chr")))
      "UCSC" else "NCBI"

    if (ann_style != frag_style) {
      if (verbose) message("[scGAS] Setting annotation seqlevels style to ", frag_style)
      GenomeInfoDb::seqlevelsStyle(annotation) <- frag_style
    } else {
      if (verbose) message("[scGAS] Annotation seqlevels style already ", frag_style, ".")
    }
  }

  if (!is.null(genome_str)) {
    if (verbose) message("[scGAS] Setting annotation genome to ", genome_str)
    GenomeInfoDb::genome(annotation) <- genome_str
  }

  return(annotation)
}
