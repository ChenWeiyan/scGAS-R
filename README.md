# scGAS <img src="man/figures/logo.png" align="right" height="139" />

**Single-Cell Gene Activation Potential Inference via a Gene Regulatory Reference Map**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`scGAS` infers **Gene Activation Potential (GAS)** from single-cell ATAC-seq
data by combining:

- A **reference map** of cCRE-Gene associations built from 167 paired ENCODE
  bulk DNase-seq / RNA-seq samples.
- **Per-gene Lasso regression models** that quantify the contribution of each
  cis-regulatory element to gene activation.
- A **Metacell + network propagation** strategy that overcomes the extreme
  sparsity of scATAC-seq data.
- A **Chromatin Potential Field** that predicts differentiation trajectories by
  linking chromatin state to RNA expression.

The method is described in:

> 

---

## Installation

```r
# Required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "Signac", "GenomicRanges", "rtracklayer",
  "EnsDb.Hsapiens.v75", "SCAVENGE"
))

install.packages(c("Seurat", "glmnet", "glmnetUtils", "igraph", "ggplot2"))

# Install scGAS from GitHub
remotes::install_github("ChenWeiyan/scGAS-R")
```

---

## Reference data

Pre-computed ENCODE reference files (~570 MB for hg19) are hosted on *xxx* and can be
downloaded in one call after installing the package:

```r
library(scGAS)

# Download all hg19 reference files into ./data/
scgas_download_reference(genome = "hg19", data_dir = "data")

# hg38 (optional):
# scgas_download_reference(genome = "hg38", data_dir = "data")
```

| File | Size | Description |
|------|------|-------------|
| `hg19_500bp_CRE.bed` | 68 MB | Reference cCRE BED (1.8 M × 500 bp bins) |
| `ENCODE_167Tissue_SigAss_pVal.05.gz` | 30 MB | cCRE–Gene association list (p < 0.05) |
| `ENCODE_167Tissue_RNA_LogNorm.gz` | 13 MB | ENCODE bulk RNA-seq, log-normalised |
| `ENCODE_167Tissue_CRE_LogNorm.gz` | 448 MB | ENCODE bulk DNase-seq, log-normalised |

---

## Quick Start

Gene annotations are loaded automatically from the matching `EnsDb` package
based on the detected genome — no manual annotation step needed.

```r
library(scGAS)

# Step 1: Preprocess scATAC-seq
# Genome and annotation are resolved automatically from the fragment header.
# out_dir / run_name are stored in obj@misc and inherited by all downstream steps.
obj <- scgas_preprocess(
  fragment_path = "fragments.tsv.gz",
  data_dir      = "data",
  out_dir       = "results",      # optional: enables checkpointing
  run_name      = "my_run",       # all output goes to results/my_run/
  n_cores       = 8
)

# Step 2: Construct Metacells
obj <- scgas_metacell(obj, n_cores = 8)

# Step 3: Train per-gene Lasso models (reference files resolved from obj@misc)
obj <- scgas_train_models(obj, n_cores = 8)

# Step 4: Compute single-cell GAS
obj <- scgas_compute(obj, n_cores = 8)

# Step 5: Chromatin Potential Field
cpf_result <- scgas_chromatin_potential(
  seurat_obj   = obj,
  point_colour = "mc_membership"
)
plot(cpf_result)
```

---

## Checkpointing & output

All pipeline functions accept three shared parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `out_dir` | `NULL` | Base output directory. Set once in `scgas_preprocess()`; all downstream functions inherit it from `obj@misc`. |
| `run_name` | `"scgas_run"` | Run label. Files are written to `out_dir/run_name/`. |
| `save_obj` | `FALSE` | Save the returned object (SeuratObject or `cpf`) to disk after each step. |

Set `out_dir` and `run_name` in `scgas_preprocess()` — they propagate forward
automatically through every subsequent call:

```r
obj <- scgas_preprocess(
  fragment_path = "fragments.tsv.gz",
  data_dir      = "data",
  out_dir       = "results",
  run_name      = "pbmc_run1",
  save_obj      = TRUE          # → results/pbmc_run1/seurat_preprocessed.rds
)

obj <- scgas_metacell(obj, save_obj = TRUE)          # → seurat_metacell.rds
obj <- scgas_train_models(obj, n_cores = 8)          # checkpoints each gene model
obj <- scgas_compute(obj, save_obj = TRUE)           # → seurat_scgas_computed.rds
cpf <- scgas_chromatin_potential(obj, save_obj = TRUE)  # → cpf.rds
```

If the run is interrupted, re-run the same commands with the same
`out_dir`/`run_name` — completed checkpoints are loaded from disk and skipped
automatically.

**Files written to `out_dir/run_name/`:**

| File | Written by | When |
|------|-----------|------|
| `seurat_preprocessed.rds` | `scgas_preprocess()` | `save_obj = TRUE` |
| `mc_membership.rds` | `scgas_metacell()` | always (checkpoint) |
| `seurat_metacell.rds` | `scgas_metacell()` | `save_obj = TRUE` |
| `atac_cl_lognorm.rds` | `scgas_train_models()` | always (checkpoint) |
| `gene_cres.rds` | `scgas_train_models()` | always (checkpoint) |
| `models/{gene}.rds` | `scgas_train_models()` | always (per-gene checkpoint) |
| `seurat_models_trained.rds` | `scgas_train_models()` | `save_obj = TRUE` |
| `scgas_mat.rds` | `scgas_compute()` | always (checkpoint) |
| `seurat_scgas_computed.rds` | `scgas_compute()` | `save_obj = TRUE` |
| `cpf.rds` | `scgas_chromatin_potential()` | `save_obj = TRUE` |

You can override `out_dir` or `run_name` in any individual call (e.g. to start
a new run from a saved checkpoint):

```r
# Load the trained object and continue with a different run name
obj <- readRDS("results/pbmc_run1/seurat_models_trained.rds")
obj <- scgas_compute(obj, out_dir = "results", run_name = "pbmc_run2", save_obj = TRUE)
```

---

## Workflow summary

Each function from Step 2 onwards takes a `SeuratObject` as its first argument
and returns the same object enriched with new results.  Step 5 is the only
exception — it returns a standalone `cpf` object.

| Step | Function | What is added |
|------|----------|---------------|
| 0 | `scgas_download_reference()` | downloads reference files from Zenodo |
| 1 | `scgas_preprocess()` | creates `SeuratObject` (ATAC assay, LSI, UMAP, KNN graph); genome and annotation auto-detected; `out_dir`/`run_name` stored in `obj@misc` |
| 2 | `scgas_metacell()` | `obj$mc_membership` + `obj@misc$mc_membership`; `mc_membership.rds` checkpoint |
| 3 | `scgas_train_models()` | `obj@misc$scgas_models`; per-gene model checkpoints in `models/` |
| 4 | `scgas_compute()` | `scGAS` assay, `scgaspca` and `scgasumap` reductions; `scgas_mat.rds` checkpoint |
| 5 | `scgas_chromatin_potential()` | returns a `cpf` object |

---

## Main functions

| Function | Description |
|----------|-------------|
| `scgas_download_reference()` | Download ENCODE reference files from Zenodo |
| `scgas_preprocess()` | Fragment → cCRE count matrix, QC, LSI, high-res clustering |
| `scgas_metacell()` | Metacell aggregation via walktrap community detection |
| `scgas_train_models()` | Per-gene Lasso model training on ENCODE bulk data |
| `scgas_compute()` | Single-cell GAS via initialisation + network propagation |
| `scgas_add_assay()` | Embed a scGAS matrix into a Seurat object |
| `scgas_chromatin_potential()` | Build and plot the Chromatin Potential Field |

---

## Citation

If you use scGAS please cite:

```
```

---

## License

MIT © 2026 Weiyan Chen
