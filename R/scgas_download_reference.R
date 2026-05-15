## ── Zenodo DOI ───────────────────────────────────────────────────────────────
## Update this single line when the Zenodo record is published.
.SCGAS_ZENODO_DOI <- NULL   # e.g. "10.5281/zenodo.1234567"


## ── File manifest ─────────────────────────────────────────────────────────────
## Lists every file hosted on Zenodo, keyed by genome assembly.
## RNA LogNorm is genome-agnostic and shared across assemblies.
.SCGAS_REFERENCE_FILES <- list(
  hg19 = c(
    "hg19_500bp_CRE.bed",
    "ENCODE_167Tissue_SigAss_pVal.05.gz",
    "ENCODE_167Tissue_RNA_LogNorm.gz",
    "ENCODE_167Tissue_CRE_LogNorm.gz"
  ),
  hg38 = c(
    "hg38_500bp_CRE.bed",
    "ENCODE_167Tissue_SigAss_hg38_pVal.05.gz",
    "ENCODE_167Tissue_RNA_LogNorm.gz",
    "ENCODE_167Tissue_CRE_hg38_LogNorm.gz"
  )
)


#' Download scGAS Reference Files from Zenodo
#'
#' Downloads the reference files needed to run the scGAS pipeline for a given
#' genome assembly and saves them into \code{data_dir}.  Already-present files
#' are skipped unless \code{overwrite = TRUE}.
#'
#' Reference files are hosted on Zenodo.  The DOI is embedded in the package
#' and will be activated when the dataset is published.  Until then, pass your
#' own \code{zenodo_doi} to use a pre-release record.
#'
#' ## Files downloaded per genome
#'
#' | File | Size (approx.) | Description |
#' |------|---------------|-------------|
#' | \code{{genome}_500bp_CRE.bed} | 70 MB | Reference cCRE BED |
#' | \code{ENCODE_167Tissue_SigAss_{genome}_pVal.05.gz} | 30 MB | cCRE–Gene associations |
#' | \code{ENCODE_167Tissue_RNA_LogNorm.gz} | 13 MB | ENCODE bulk RNA-seq |
#' | \code{ENCODE_167Tissue_CRE_{genome}_LogNorm.gz} | 450 MB | ENCODE bulk DNase-seq |
#'
#' @param genome Genome assembly to download.  Currently \code{"hg19"} and
#'   \code{"hg38"} are supported.  Default \code{"hg19"}.
#' @param data_dir Directory to save files into.  Created if it does not exist.
#'   Default \code{"data"}.
#' @param zenodo_doi Zenodo DOI string, e.g. \code{"10.5281/zenodo.1234567"}.
#'   \code{NULL} (default) uses the DOI embedded in the package.
#' @param overwrite Logical.  Re-download and overwrite files that already
#'   exist.  Default \code{FALSE}.
#' @param timeout Download timeout in seconds per file.  Default \code{3600}
#'   (1 hour), needed for the large DNase matrix.
#' @param verbose Logical.  Default \code{TRUE}.
#'
#' @return Invisibly returns a named character vector of downloaded file paths.
#'
#' @examples
#' \dontrun{
#' ## Download hg19 reference files into ./data/
#' scgas_download_reference(genome = "hg19", data_dir = "data")
#'
#' ## Download hg38 files once the DOI is live
#' scgas_download_reference(genome = "hg38", data_dir = "data")
#'
#' ## Use a pre-release or custom Zenodo record
#' scgas_download_reference(genome = "hg19", data_dir = "data",
#'                          zenodo_doi = "10.5281/zenodo.1234567")
#' }
#'
#' @export
scgas_download_reference <- function(genome      = "hg19",
                                     data_dir    = "data",
                                     zenodo_doi  = NULL,
                                     overwrite   = FALSE,
                                     timeout     = 3600L,
                                     verbose     = TRUE) {

  ## ── Resolve DOI ────────────────────────────────────────────────────────────
  doi <- if (!is.null(zenodo_doi)) zenodo_doi else .SCGAS_ZENODO_DOI

  if (is.null(doi)) {
    stop(
      "[scGAS] Zenodo DOI not yet set.\n",
      "  The reference dataset has not been published yet.\n",
      "  Check https://github.com/ChenWeiyan/scGAS-R for updates, or\n",
      "  pass zenodo_doi = '10.5281/zenodo.XXXXXXX' if you have a pre-release link."
    )
  }

  ## ── Resolve genome ─────────────────────────────────────────────────────────
  supported <- names(.SCGAS_REFERENCE_FILES)
  if (!genome %in% supported)
    stop("[scGAS] Unsupported genome '", genome, "'. ",
         "Available: ", paste(supported, collapse = ", "))

  files <- .SCGAS_REFERENCE_FILES[[genome]]

  ## ── Build base URL from DOI ────────────────────────────────────────────────
  ## Zenodo DOI format: "10.5281/zenodo.XXXXXXX"
  ## Download URL:      "https://zenodo.org/record/XXXXXXX/files/filename"
  record_id <- sub("^10\\.5281/zenodo\\.", "", doi)
  base_url  <- sprintf("https://zenodo.org/record/%s/files", record_id)

  ## ── Create data_dir if needed ──────────────────────────────────────────────
  if (!dir.exists(data_dir)) {
    if (verbose) message("[scGAS] Creating directory: ", data_dir)
    dir.create(data_dir, recursive = TRUE)
  }
  data_dir <- normalizePath(data_dir)

  ## ── Download each file ─────────────────────────────────────────────────────
  dest_paths <- stats::setNames(file.path(data_dir, files), files)

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = timeout)

  for (fname in files) {
    dest <- dest_paths[[fname]]

    if (file.exists(dest) && !overwrite) {
      if (verbose) message("[scGAS] Skipping (already exists): ", fname)
      next
    }

    url <- paste0(base_url, "/", fname, "?download=1")
    if (verbose) message("[scGAS] Downloading: ", fname,
                         "  (", .format_size(url), ")")

    tryCatch(
      utils::download.file(url, destfile = dest, mode = "wb", quiet = !verbose),
      error = function(e) {
        if (file.exists(dest)) unlink(dest)
        stop("[scGAS] Failed to download '", fname, "':\n  ", conditionMessage(e))
      }
    )

    if (verbose) message("[scGAS] Saved: ", dest)
  }

  if (verbose) message("[scGAS] All reference files ready in: ", data_dir)
  invisible(dest_paths)
}


## ── Helper: human-readable file size (best-effort from Content-Length) ───────
#' @keywords internal
.format_size <- function(url) {
  tryCatch({
    h <- curlGetHeaders(url, redirect = TRUE, verify = FALSE)
    cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
    if (length(cl) == 0) return("size unknown")
    bytes <- as.numeric(trimws(sub(".*:", "", cl[1])))
    if (bytes >= 1e9) sprintf("%.1f GB", bytes / 1e9)
    else if (bytes >= 1e6) sprintf("%.0f MB", bytes / 1e6)
    else if (bytes >= 1e3) sprintf("%.0f KB", bytes / 1e3)
    else sprintf("%d B", as.integer(bytes))
  }, error = function(e) "size unknown")
}
