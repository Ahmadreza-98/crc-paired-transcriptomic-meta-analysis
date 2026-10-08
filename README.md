# CRC paired-cohort transcriptomic meta-analysis

This repository contains the R analysis pipeline for:

**Reproducible Downregulation and Network Convergence of AQP8, GUCA2A, and MS4A12 in Colorectal Cancer: A Paired-Cohort Transcriptomic Meta-analysis**

**Authors:** Seyed Ahmadreza Siadat, Eghbal Mansoori, Navid Mogharrab, Mostafa Saadat  
**Maintainer:** Seyed Ahmadreza Siadat ([ORCID 0009-0004-1439-6306](https://orcid.org/0009-0004-1439-6306))

The pipeline analyzes nine paired colorectal cancer GEO cohorts, performs gene-wise random-effects meta-analysis, enrichment and network analyses, and evaluates AQP8, GUCA2A, and MS4A12 in TCGA-COAD/READ.

## Repository contents

```text
crc-paired-transcriptomic-meta-analysis/
├── analysis.R
├── README.md
├── DATA_PROVENANCE.md
├── CITATION.cff
├── LICENSE
├── crc-paired-transcriptomic-meta-analysis.Rproj
├── Data/
│   └── README.md
├── Cache/
│   └── README.md
├── Results_manuscript/
│   └── README.md
└── environment/
    ├── README.md
    └── sessionInfo_tested.txt
```

Large input files and complete generated results are distributed separately as GitHub Release assets rather than stored in Git history.

## Analysis overview

`analysis.R` performs the full analysis reported in the manuscript:

- paired limma analysis within each GEO cohort;
- REML random-effects meta-analysis across cohorts;
- leave-one-study-out sensitivity analysis for the focal genes;
- GO, KEGG, Hallmark, cytoband, and regulatory-target enrichment;
- cross-platform integration, ComBat adjustment, and signed WGCNA;
- STRING v12.0 PPI reconstruction and centrality analysis;
- TCGA-COAD/READ expression validation using STAR counts and TPM;
- adjusted Cox models, expression-by-stage interaction tests, and Kaplan-Meier plots.

The focal genes are **AQP8, GUCA2A, and MS4A12**.

## Data

The GEO analysis uses these nine series:

`GSE184093`, `GSE156355`, `GSE89076`, `GSE84984`, `GSE75970`, `GSE74602`, `GSE22598`, `GSE25070`, and `GSE21510`.

For the manuscript reproduction workflow, download the GEO/GPL input archive from the repository's **GitHub Releases** page and place the extracted files directly in `Data/`. The required filenames are listed in [`Data/README.md`](Data/README.md).

TCGA-COAD/READ, STRING, and other online resources are retrieved by the pipeline. Their roles are summarized in [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md).

## Requirements

The study analysis was run with **R 4.4.1**. The exact R session information from the validated final run is provided in [`environment/sessionInfo_tested.txt`](environment/sessionInfo_tested.txt).

The script checks for required packages at startup and does not install or update packages automatically.

Main CRAN packages include `dplyr`, `tidyr`, `purrr`, `readr`, `tibble`, `stringr`, `data.table`, `ggplot2`, `ggrepel`, `patchwork`, `scales`, `metafor`, `WGCNA`, `pheatmap`, `igraph`, `ggraph`, `png`, `survival`, `survminer`, `msigdbr`, `matrixStats`, `httr2`, and `jsonlite`.

Main Bioconductor packages include `GEOquery`, `Biobase`, `limma`, `sva`, `AnnotationDbi`, `org.Hs.eg.db`, `clusterProfiler`, `TCGAbiolinks`, `SummarizedExperiment`, `DESeq2`, and `STRINGdb`.

## Running the analysis

1. Extract or clone the repository.
2. Add the GEO/GPL input files to `Data/`.
3. Open `crc-paired-transcriptomic-meta-analysis.Rproj` in RStudio, or start R from the repository root.
4. Run:

```r
source("analysis.R")
```

The script writes tables and figures under `Results_manuscript/` and records the R session in `Results_manuscript/sessionInfo.txt`.

TCGA RNA-seq files are cached locally. On Windows, the script may use `~/crc_tcga_cache` when a shorter path is needed for `TCGAbiolinks` downloads.

For the survival analysis, case-level `sex_at_birth` is retrieved from the GDC Cases API and used as the sex covariate. If this required query cannot be completed after retries, the pipeline stops rather than substituting the incomplete sex/gender field embedded in the RNA-seq object.

## Outputs

A complete run produces manuscript and supplementary figures, machine-readable tables, QC files, and `sessionInfo.txt` under `Results_manuscript/`.

The complete output archive from the validated manuscript run is available from the repository's **GitHub Releases** page. The tested R session information retained in `environment/sessionInfo_tested.txt` corresponds to that validated final run.

## Citation and license

Repository citation metadata are provided in `CITATION.cff`.

The code and repository documentation are released under the [MIT License](LICENSE). Third-party datasets and resources remain subject to their own terms of use.
