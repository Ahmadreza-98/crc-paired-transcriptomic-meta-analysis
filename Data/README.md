# Data

Place the GEO Series Matrix files and GPL annotation tables used by the analysis in this directory.

The easiest way to reproduce the manuscript workflow is to download the input-data archive from the repository's **GitHub Releases** page and extract its contents directly into `Data/`.

## Required GEO files

```text
GSE184093.gz
GSE156355.gz
GSE89076.gz
GSE84984.gz
GSE75970.gz
GSE74602.gz
GSE22598.gz
GSE25070.gz
GSE21510.gz
```

The standard GEO filename form `<GSE>_series_matrix.txt.gz` is also accepted by `analysis.R`.

## Required GPL files

```text
GPL20115.txt
GPL21185.txt
GPL16699.txt
GPL17586.txt
GPL14550.txt
GPL6104.txt
GPL570.txt
GPL6883.txt
```

`GPL570` is shared by GSE22598 and GSE21510.

A typical directory therefore looks like:

```text
Data/
├── README.md
├── GSE184093.gz
├── GSE156355.gz
├── GSE89076.gz
├── GSE84984.gz
├── GSE75970.gz
├── GSE74602.gz
├── GSE22598.gz
├── GSE25070.gz
├── GSE21510.gz
├── GPL20115.txt
├── GPL21185.txt
├── GPL16699.txt
├── GPL17586.txt
├── GPL14550.txt
├── GPL6104.txt
├── GPL570.txt
└── GPL6883.txt
```

The same public files can also be obtained directly from NCBI GEO. Local copies are recommended for manuscript reproduction so the GEO/GPL inputs used for the run are explicit.

`Data/PPI_centrality_export.tsv` is optional and is only used when a retained Cytoscape centrality export is supplied.

TCGA/GDC and STRING data are retrieved by `analysis.R` and are not stored in this directory.
