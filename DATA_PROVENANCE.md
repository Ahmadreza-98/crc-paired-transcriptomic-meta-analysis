# Data provenance

This file summarizes the external data sources used by `analysis.R`.

## GEO microarray data

The paired meta-analysis uses nine processed GEO Series Matrix datasets and their GEO platform annotation tables:

| GEO series | Platform | Paired patients |
| --- | --- | ---: |
| GSE184093 | GPL20115 | 9 |
| GSE156355 | GPL21185 | 6 |
| GSE89076 | GPL16699 | 34 |
| GSE84984 | GPL17586 | 5 |
| GSE75970 | GPL14550 | 3 |
| GSE74602 | GPL6104 | 30 |
| GSE22598 | GPL570 | 16 |
| GSE25070 | GPL6883 | 22 |
| GSE21510 | GPL570 | 20 |

The final GEO dataset contains **145 paired patients (290 samples)**. The predefined sample subsets, exclusions, technical-replicate handling, transformation rule, and platform-to-Entrez mapping are encoded in `analysis.R`.

For WGCNA, six predefined patient pairs are removed as complete pairs before network construction, leaving 139 paired patients (278 samples). Their global patient indices are `28, 38, 44, 48, 90, 103`.

The exact GEO/GPL input files used for the manuscript workflow are distributed as a GitHub Release asset and can be placed directly in `Data/`.

## TCGA / GDC

TCGA-COAD and TCGA-READ RNA-seq data are retrieved from the NCI Genomic Data Commons using `TCGAbiolinks`.

- STAR counts are used for DESeq2 expression validation.
- TPM values are used for focal-gene expression plots and survival analyses.
- Case-level `sex_at_birth` is retrieved from the GDC Cases API and used as the sex covariate in the survival models.

The sex/gender field embedded in the RNA-seq `SummarizedExperiment` is used only for comparison and is not used as a fallback for the survival models. If the required GDC `sex_at_birth` query fails after retries, the pipeline stops.

TCGA/GDC is a live resource. If upstream records change, a future rerun can legitimately differ from an earlier run.

## STRING

The PPI analysis uses **STRING v12.0** with a combined-score threshold of **400**. STRING files are retrieved through `STRINGdb` and cached locally.

Five MCODE module memberships retained from the earlier Cytoscape workflow are encoded in `analysis.R`; MCODE clustering is not rerun in the final R pipeline. The centrality metrics used in the manuscript are calculated in R. If available, `Data/PPI_centrality_export.tsv` can be supplied as an optional retained Cytoscape export.

## Functional enrichment

GO and KEGG analyses use `clusterProfiler`. Hallmark, cytoband, and regulatory-target gene sets are obtained through `msigdbr` and associated MSigDB-compatible resources.

Because these databases are updated over time, enrichment results from a future rerun can differ if the upstream resource versions change.

## Environment record

A successful run writes `Results_manuscript/sessionInfo.txt`. The manuscript reports **R 4.4.1** for the study analysis.
