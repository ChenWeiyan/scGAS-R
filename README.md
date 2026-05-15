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

```r
library(scGAS)
library(EnsDb.Hsapiens.v75)
library(Signac)

# Gene annotations — genome style will be set automatically by scgas_preprocess()
# Choose bsed on your data
annotation <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v75)

# Step 1: Preprocess scATAC-seq (genome auto-detected from fragment header)
obj <- scgas_preprocess(
  fragment_path     = "fragments.tsv.gz",
  data_dir          = "data",
  annotation        = annotation,
  n_cores           = 8
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

## Workflow summary

Each function from Step 2 onwards takes a `SeuratObject` as its first argument
and returns the same object enriched with new results.  Step 5 is the only
exception — it returns a standalone `cpf` object.

| Step | Function | What is added |
|------|----------|---------------|
| 0 | `scgas_download_reference()` | downloads reference files from Zenodo |
| 1 | `scgas_preprocess()` | creates `SeuratObject` (ATAC assay, LSI, UMAP, KNN graph) |
| 2 | `scgas_metacell()` | `obj$mc_membership` + `obj@misc$mc_membership` |
| 3 | `scgas_train_models()` | `obj@misc$scgas_models` |
| 4 | `scgas_compute()` | `scGAS` assay, `scgaspca` and `scgasumap` reductions |
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
