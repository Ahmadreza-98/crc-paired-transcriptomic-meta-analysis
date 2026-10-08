# Cache

`analysis.R` uses this directory for downloaded files that can be reused between runs.

```text
Cache/
├── GEO/
├── TCGA/
└── STRING/
```

- `GEO/` stores Series Matrix files obtained by the network fallback.
- `TCGA/` stores GDC/TCGA downloads and prepared objects.
- `STRING/` stores files used by `STRINGdb`.

Cache contents are ignored by Git and can be regenerated.

On Windows, `TCGAbiolinks` can create very long file paths. When needed, the script uses the shorter user directory `~/crc_tcga_cache` instead of `Cache/TCGA/`. The active path is printed at startup. A custom path can be set with `TCGA_CACHE_OVERRIDE` near the start of `analysis.R`.
