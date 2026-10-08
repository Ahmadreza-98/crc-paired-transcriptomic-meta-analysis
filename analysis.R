###############################################################################
# Colorectal cancer integrative transcriptomic analysis
# Final reproducible analysis pipeline for the manuscript
#
# The pipeline includes paired GEO differential expression, random-effects
# meta-analysis, pathway enrichment, WGCNA, STRING/PPI analysis, TCGA
# expression validation, and TPM-based survival analysis.
#
# Focal genes: AQP8, GUCA2A, MS4A12
###############################################################################

options(stringsAsFactors = FALSE)
options(timeout = max(600, getOption("timeout")))
options(download.file.method = "libcurl")
set.seed(1234)

###############################################################################
# 0. PROJECT SETTINGS
###############################################################################

PROJECT_ROOT_OVERRIDE <- NULL

normalize_project_candidate <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(NULL)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

is_project_root <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  file.exists(file.path(path, "analysis.R")) &&
    dir.exists(file.path(path, "Data"))
}

ancestor_directories <- function(path, max_levels = 12L) {
  path <- normalize_project_candidate(path)
  if (is.null(path)) return(character(0))

  # If a file path was supplied, start from its containing directory.
  if (file.exists(path) && !dir.exists(path)) {
    path <- dirname(path)
  }

  out <- character(0)
  current <- path

  for (i in seq_len(max_levels)) {
    out <- c(out, current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }

  unique(out)
}

detect_source_directory <- function() {
  # source("analysis.R") keeps the sourced filename in one of the active frames.
  frame_files <- unlist(
    lapply(
      sys.frames(),
      function(fr) {
        x <- tryCatch(fr$ofile, error = function(e) NULL)
        if (is.null(x) || length(x) != 1L || !nzchar(x)) {
          character(0)
        } else {
          as.character(x)
        }
      }
    ),
    use.names = FALSE
  )

  if (length(frame_files) == 0L) return(NULL)

  frame_files <- vapply(
    frame_files,
    function(x) normalizePath(x, winslash = "/", mustWork = FALSE),
    character(1)
  )

  analysis_file <- frame_files[basename(frame_files) == "analysis.R"]
  chosen_file <- if (length(analysis_file) > 0L) {
    tail(analysis_file, 1L)
  } else {
    tail(frame_files, 1L)
  }

  dirname(chosen_file)
}

detect_rscript_directory <- function() {
  # Rscript analysis.R exposes the script path through the --file argument.
  command_line <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_line, value = TRUE)

  if (length(file_arg) == 0L) return(NULL)

  script_file <- sub("^--file=", "", file_arg[1L])
  script_file <- normalizePath(
    script_file,
    winslash = "/",
    mustWork = FALSE
  )
  dirname(script_file)
}

detect_rstudio_directories <- function() {
  # When the whole script is sent to the Console with RStudio's Run command,
  # source() metadata are not available.  In that situation, use the active
  # RStudio document/project path when the optional rstudioapi package exists.
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    return(character(0))
  }

  available <- tryCatch(
    rstudioapi::isAvailable(),
    error = function(e) FALSE
  )
  if (!isTRUE(available)) return(character(0))

  out <- character(0)

  active_file <- tryCatch(
    rstudioapi::getActiveDocumentContext()$path,
    error = function(e) ""
  )
  if (length(active_file) == 1L && nzchar(active_file)) {
    out <- c(out, dirname(active_file))
  }

  active_project <- tryCatch(
    rstudioapi::getActiveProject(),
    error = function(e) ""
  )
  if (length(active_project) == 1L && nzchar(active_project)) {
    out <- c(out, active_project)
  }

  unique(out)
}

find_project_root <- function(override = NULL) {
  raw_candidates <- c(
    override,
    detect_source_directory(),
    detect_rscript_directory(),
    detect_rstudio_directories(),
    getwd()
  )

  raw_candidates <- unique(
    raw_candidates[
      !is.na(raw_candidates) & nzchar(raw_candidates)
    ]
  )

  candidates <- unique(
    unlist(
      lapply(raw_candidates, ancestor_directories),
      use.names = FALSE
    )
  )

  hits <- candidates[
    vapply(candidates, is_project_root, logical(1))
  ]

  if (length(hits) > 0L) {
    return(hits[1L])
  }

  # One additional safe fallback: if the current working directory is the
  # parent of the repository, accept a unique direct child that looks like it.
  child_dirs <- tryCatch(
    list.dirs(getwd(), recursive = FALSE, full.names = TRUE),
    error = function(e) character(0)
  )
  child_hits <- child_dirs[
    vapply(child_dirs, is_project_root, logical(1))
  ]
  if (length(child_hits) == 1L) {
    return(normalizePath(child_hits[1L], winslash = "/", mustWork = FALSE))
  }

  # If code was pasted/run line-by-line and no file context is available, an
  # interactive session cannot infer the repository location reliably.  Give
  # the user a final chance to select it instead of failing immediately.
  if (interactive()) {
    selected_dir <- NULL

    if (.Platform$OS.type == "windows") {
      selected_dir <- tryCatch(
        utils::choose.dir(
          default = getwd(),
          caption = "Select the crc-paired-transcriptomic-meta-analysis folder"
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(selected_dir) && length(selected_dir) == 1L &&
        !is.na(selected_dir) && nzchar(selected_dir) &&
        is_project_root(selected_dir)) {
      return(normalizePath(selected_dir, winslash = "/", mustWork = FALSE))
    }
  }

  stop(
    "Project directory not found.\n",
    "Open the included .Rproj and run analysis.R, use source('analysis.R') from ",
    "the repository, or set PROJECT_ROOT_OVERRIDE to the repository path.\n",
    "The project root must contain analysis.R and the Data directory."
  )
}

PROJECT_ROOT <- find_project_root(PROJECT_ROOT_OVERRIDE)
message("Project root: ", PROJECT_ROOT)

DIR_DATA          <- file.path(PROJECT_ROOT, "Data")
DIR_RESULTS       <- file.path(PROJECT_ROOT, "Results_manuscript")
DIR_FIG           <- file.path(DIR_RESULTS, "Figures")
DIR_TABLE         <- file.path(DIR_RESULTS, "Tables")
DIR_MANUSCRIPT    <- file.path(DIR_FIG, "Manuscript")
DIR_SUPPLEMENTARY <- file.path(DIR_FIG, "Supplementary")
DIR_CACHE         <- file.path(PROJECT_ROOT, "Cache")
DIR_GEO           <- file.path(DIR_CACHE, "GEO")
DIR_STRING        <- file.path(DIR_CACHE, "STRING")

# TCGAbiolinks expands downloaded GDC files into deeply nested directories.
# On Windows, a long repository path can make those generated file paths too
# long for some file operations. Keep the normal project cache when the path is
# short; otherwise use a short persistent per-user TCGA cache directory.
TCGA_CACHE_OVERRIDE <- NULL
TCGA_PROJECT_CACHE <- file.path(DIR_CACHE, "TCGA")

if (!is.null(TCGA_CACHE_OVERRIDE) &&
    length(TCGA_CACHE_OVERRIDE) == 1L &&
    !is.na(TCGA_CACHE_OVERRIDE) &&
    nzchar(TCGA_CACHE_OVERRIDE)) {
  DIR_TCGA <- normalizePath(
    TCGA_CACHE_OVERRIDE,
    winslash = "/",
    mustWork = FALSE
  )
} else if (
  .Platform$OS.type == "windows" &&
    nchar(normalizePath(TCGA_PROJECT_CACHE, winslash = "/", mustWork = FALSE)) > 60L
) {
  DIR_TCGA <- file.path(path.expand("~"), "crc_tcga_cache")
} else {
  DIR_TCGA <- TCGA_PROJECT_CACHE
}

message("TCGA cache: ", DIR_TCGA)

invisible(
  lapply(
    c(
      DIR_RESULTS,
      DIR_FIG,
      DIR_TABLE,
      DIR_MANUSCRIPT,
      DIR_SUPPLEMENTARY,
      DIR_CACHE,
      DIR_GEO,
      DIR_TCGA,
      DIR_STRING
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)

META_FDR_CUTOFF   <- 0.05
META_LOGFC_CUTOFF <- 1.5
META_MIN_STUDIES  <- 3
META_METHOD        <- "REML"

WGCNA_SFT_R2_TARGET     <- 0.85
WGCNA_MIN_MODULE_SIZE   <- 50
WGCNA_MERGE_CUT         <- 0.20
WGCNA_MM_CUTOFF         <- 0.85
WGCNA_GS_CUTOFF         <- 0.85

WGCNA_OUTLIER_GLOBAL_PATIENTS <- c(28, 38, 44, 48, 90, 103)

FOCAL_GENES <- c("AQP8", "GUCA2A", "MS4A12")

STRING_VERSION <- "12.0"
STRING_SCORE_THRESHOLD <- 400

PLOT_COLORS <- c(
  down = "#2C7BB6",
  up = "#D7191C",
  normal = "#2C7BB6",
  tumor = "#D7191C",
  early = "#2C7BB6",
  advanced = "#D7191C",
  low = "#2C7BB6",
  high = "#D7191C",
  accent = "#6A51A3",
  neutral = "#BDBDBD"
)

GENE_COLORS <- c(
  AQP8 = "#2C7BB6",
  GUCA2A = "#1B9E77",
  MS4A12 = "#D95F02"
)

###############################################################################
# 1. PACKAGES
###############################################################################

required_packages <- c(
  # CRAN
  "dplyr", "tidyr", "purrr", "readr", "tibble", "stringr", "data.table",
  "ggplot2", "ggrepel", "patchwork", "scales", "metafor",
  "WGCNA", "pheatmap", "igraph", "ggraph", "png",
  "survival", "survminer", "msigdbr", "matrixStats", "httr2", "jsonlite",
  # Bioconductor
  "GEOquery", "Biobase", "limma", "sva", "AnnotationDbi",
  "org.Hs.eg.db", "clusterProfiler", "TCGAbiolinks",
  "SummarizedExperiment", "DESeq2", "STRINGdb"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "), "\n\n",
    "Install CRAN packages with install.packages(...).\n",
    "Install Bioconductor packages with BiocManager::install(...).\n",
    "Then rerun the script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

WGCNA::allowWGCNAThreads()

###############################################################################
# 2. PLOT STYLE AND GENERAL UTILITIES
###############################################################################

theme_set(
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = 13)
    )
)

save_plot_to <- function(
    plot,
    directory,
    filename,
    width = 7,
    height = 5,
    dpi = 400
) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    file.path(directory, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
  
  ggsave(
    file.path(directory, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
}

save_plot <- function(plot, filename, width = 7, height = 5, dpi = 400) {
  save_plot_to(plot, DIR_FIG, filename, width, height, dpi)
}

save_manuscript_plot <- function(plot, filename, width = 7, height = 5, dpi = 600) {
  save_plot_to(plot, DIR_MANUSCRIPT, filename, width, height, dpi)
}

save_supplementary_plot <- function(plot, filename, width = 7, height = 5, dpi = 500) {
  save_plot_to(plot, DIR_SUPPLEMENTARY, filename, width, height, dpi)
}

write_tsv_safe <- function(x, filename) {
  readr::write_tsv(
    as.data.frame(x),
    file.path(DIR_TABLE, filename),
    na = "NA"
  )
}

write_matrix_tsv <- function(x, filename) {
  out <- as.data.frame(x) |>
    tibble::rownames_to_column("gene_id")
  
  readr::write_tsv(
    out,
    file.path(DIR_TABLE, filename),
    na = "NA"
  )
}


retry_remote_call <- function(
    label,
    fun,
    attempts = 4L,
    wait_seconds = c(10, 30, 60)
) {
  attempts <- max(1L, as.integer(attempts))
  last_error <- NULL
  
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(
      fun(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    
    if (!is.null(result)) {
      return(result)
    }
    
    if (attempt < attempts) {
      wait_now <- wait_seconds[
        min(attempt, length(wait_seconds))
      ]
      
      message(
        label,
        " failed on attempt ",
        attempt,
        "/",
        attempts,
        ": ",
        conditionMessage(last_error),
        "\nRetrying in ",
        wait_now,
        " seconds ..."
      )
      
      Sys.sleep(wait_now)
    }
  }
  
  stop(
    label,
    " failed after ",
    attempts,
    " attempts.\nLast error: ",
    conditionMessage(last_error),
    "\nThe analysis requires this external service for a complete manuscript run.",
    call. = FALSE
  )
}

format_p <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "NA",
    format.pval(x, digits = digits, eps = 0.001)
  )
}

minmax01 <- function(x) {
  x <- as.numeric(x)
  if (all(!is.finite(x))) return(rep(0, length(x)))
  rng <- range(x[is.finite(x)], na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(1, length(x)))
  (x - rng[1]) / diff(rng)
}

ratio_to_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2) return(NA_real_)
    as.numeric(z[1]) / as.numeric(z[2])
  }, numeric(1))
}

needs_log2 <- function(mat) {
  x <- as.numeric(mat)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(FALSE)
  q <- stats::quantile(x, c(0.01, 0.25, 0.50, 0.75, 0.99), na.rm = TRUE)
  isTRUE(q[5] > 100 || (q[5] - q[1] > 50 && q[2] > 0))
}

log2_if_needed <- function(mat) {
  if (needs_log2(mat)) log2(mat + 1) else mat
}

# The uploaded GPL files were inspected directly.  For the final nine cohorts,
# every platform supplies an Entrez/LocusLink identifier explicitly (GPL17586 via
# the fifth field of gene_assignment).  We therefore use the platform-provided
# Entrez identifier as the cross-platform key, exactly as in the original analysis.
# This is preferable here to remapping through current gene symbols, because symbol
# aliases/database updates can change the historical probe-to-gene universe.
first_entrez <- function(x) {
  # Reproduce the identifier rule used in the original analysis exactly:
  # keep the FIRST whitespace-delimited token from the dedicated Entrez/LocusLink
  # annotation field and retain it only when it is a pure integer.
  #
  # This is intentionally different from searching for the first number anywhere
  # in the string.  For GPL17586, a first assignment of "---" followed by a later
  # assignment must remain unmapped, matching the historical strsplit2(...)[,5]
  # followed by sub(' .*','', ...).  This exact rule reproduces the validated
  # 15,890-gene intersection across the uploaded nine GEO cohorts.
  x <- as.character(x)
  x[is.na(x)] <- ""
  token <- sub(" .*", "", x)
  token <- trimws(token)
  token[token %in% c("", "---", "NA", "na")] <- NA_character_
  token[!is.na(token) & !grepl("^[0-9]+$", token)] <- NA_character_
  token
}

first_symbol <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "---", "NA", "na")] <- NA_character_
  x <- sub("\\s*///.*$", "", x)
  x <- sub("\\s*//.*$", "", x)
  x <- sub(";.*$", "", x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

find_annotation_col <- function(df, candidates) {
  normalize_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  cn <- colnames(df)
  hit <- match(normalize_name(candidates), normalize_name(cn), nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) == 0) return(NA_character_)
  cn[hit[1]]
}

aggregate_probes_to_entrez <- function(expr, entrez) {
  ok <- !is.na(entrez) & grepl("^[0-9]+$", entrez)
  expr <- expr[ok, , drop = FALSE]
  entrez <- entrez[ok]
  if (nrow(expr) == 0) stop("No probes could be mapped to Entrez IDs.")
  
  sums <- rowsum(expr, group = entrez, reorder = FALSE, na.rm = TRUE)
  n_per_gene <- as.numeric(table(factor(entrez, levels = rownames(sums))))
  out <- sweep(sums, 1, n_per_gene, "/")
  out
}

make_condition <- function(pattern, n_pairs) {
  switch(
    pattern,
    alt_tumor_normal = rep(c("tumor", "normal"), n_pairs),
    alt_normal_tumor = rep(c("normal", "tumor"), n_pairs),
    block_normal_tumor = c(rep("normal", n_pairs), rep("tumor", n_pairs)),
    block_tumor_normal = c(rep("tumor", n_pairs), rep("normal", n_pairs)),
    stop("Unknown pairing pattern: ", pattern)
  )
}

make_patient_id <- function(pattern, n_pairs) {
  if (grepl("^alt_", pattern)) {
    rep(seq_len(n_pairs), each = 2)
  } else {
    rep(seq_len(n_pairs), 2)
  }
}

check_pairs <- function(meta, label) {
  chk <- meta |>
    count(patient_local, condition) |>
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)
  if (!all(c("normal", "tumor") %in% colnames(chk))) {
    stop(label, ": pairing check failed; both normal and tumor columns are required.")
  }
  bad <- chk |> filter(normal != 1 | tumor != 1)
  if (nrow(bad) > 0) {
    stop(label, ": at least one patient does not have exactly one normal and one tumor sample.")
  }
  invisible(TRUE)
}

###############################################################################
# 3. GEO COHORT CONFIGURATION
#    This reproduces the FINAL NINE cohorts and the original sample exclusions.
###############################################################################

GSE89076_SAMPLES <- c(
  "GSM2358437", "GSM2358438", "GSM2358439", "GSM2358440", "GSM2358441", "GSM2358442",
  "GSM2358443", "GSM2358444", "GSM2358445", "GSM2358446", "GSM2358447", "GSM2358508",
  "GSM2358448", "GSM2358449", "GSM2358450", "GSM2358451", "GSM2358509", "GSM2358452",
  "GSM2358453", "GSM2358454", "GSM2358510", "GSM2358455", "GSM2358456", "GSM2358457",
  "GSM2358458", "GSM2358459", "GSM2358504", "GSM2358460", "GSM2358505", "GSM2358461",
  "GSM2358462", "GSM2358463", "GSM2358464", "GSM2358465", "GSM2358466", "GSM2358467",
  "GSM2358468", "GSM2358469", "GSM2358514", "GSM2358470", "GSM2358501", "GSM2358473",
  "GSM2358502", "GSM2358474", "GSM2358475", "GSM2358476", "GSM2358477", "GSM2358478",
  "GSM2358479", "GSM2358480", "GSM2358515", "GSM2358481", "GSM2358482", "GSM2358483",
  "GSM2358484", "GSM2358485", "GSM2358486", "GSM2358487", "GSM2358488", "GSM2358489",
  "GSM2358506", "GSM2358490", "GSM2358491", "GSM2358492", "GSM2358507", "GSM2358493",
  "GSM2358494", "GSM2358495", "GSM2358496", "GSM2358497", "GSM2358498", "GSM2358499",
  "GSM2358516", "GSM2358500", "GSM2358503", "GSM2358511", "GSM2358512", "GSM2358513"
)

GSE22598_SAMPLES <- c(
  "GSM452629", "GSM452630", "GSM452631", "GSM452632", "GSM452633", "GSM452634",
  "GSM452635", "GSM452636", "GSM452637", "GSM452638", "GSM452639", "GSM452640",
  "GSM452641", "GSM452642", "GSM452643", "GSM452644", "GSM452645", "GSM452646",
  "GSM452647", "GSM452648", "GSM452649", "GSM452650", "GSM452651", "GSM452652",
  "GSM452653", "GSM452654", "GSM452655", "GSM452656", "GSM452657", "GSM452658",
  "GSM452659", "GSM452660", "GSM452661", "GSM452662"
)

GSE21510_SAMPLES <- c(
  "GSM549099", "GSM537336", "GSM549100", "GSM537338", "GSM549101", "GSM537339",
  "GSM549102", "GSM537341", "GSM549103", "GSM537343", "GSM549104", "GSM537345",
  "GSM549105", "GSM537351", "GSM549106", "GSM537352", "GSM549107", "GSM537353",
  "GSM549108", "GSM537355", "GSM549109", "GSM537356", "GSM549110", "GSM537383",
  "GSM549111", "GSM549142", "GSM549112", "GSM537384", "GSM549113", "GSM537386",
  "GSM549114", "GSM537388", "GSM549115", "GSM537389", "GSM549116", "GSM537390",
  "GSM549117", "GSM537391", "GSM549119", "GSM537392", "GSM549120", "GSM537393",
  "GSM549122", "GSM537395", "GSM549123", "GSM537397"
)

geo_cfg <- list(
  GSE184093 = list(
    gse = "GSE184093", gpl = "GPL20115", pattern = "alt_tumor_normal",
    n_pairs = 9, sample_subset = NULL, drop_indices = integer(0),
    technical_replicates = list(), expected_final_pairs = 9
  ),
  GSE156355 = list(
    gse = "GSE156355", gpl = "GPL21185", pattern = "block_normal_tumor",
    n_pairs = 6, sample_subset = NULL, drop_indices = integer(0),
    technical_replicates = list(), expected_final_pairs = 6
  ),
  GSE89076 = list(
    gse = "GSE89076", gpl = "GPL16699", pattern = "alt_normal_tumor",
    n_pairs = 39, sample_subset = GSE89076_SAMPLES,
    drop_indices = c(27, 28, 49, 50, 59, 60, 65, 66, 75, 76),
    technical_replicates = list(), expected_final_pairs = 34
  ),
  GSE84984 = list(
    gse = "GSE84984", gpl = "GPL17586", pattern = "block_tumor_normal",
    n_pairs = 6, sample_subset = NULL, drop_indices = c(6, 12),
    technical_replicates = list(
      c("GSM2255446", "GSM2255447"),
      c("GSM2255449", "GSM2255450"),
      c("GSM2255452", "GSM2255453")
    ), expected_final_pairs = 5
  ),
  GSE75970 = list(
    gse = "GSE75970", gpl = "GPL14550", pattern = "alt_tumor_normal",
    n_pairs = 4, sample_subset = NULL, drop_indices = c(5, 6),
    technical_replicates = list(), expected_final_pairs = 3
  ),
  GSE74602 = list(
    gse = "GSE74602", gpl = "GPL6104", pattern = "alt_tumor_normal",
    n_pairs = 30, sample_subset = NULL, drop_indices = integer(0),
    technical_replicates = list(), expected_final_pairs = 30
  ),
  GSE22598 = list(
    gse = "GSE22598", gpl = "GPL570", pattern = "block_normal_tumor",
    n_pairs = 17, sample_subset = GSE22598_SAMPLES, drop_indices = c(11, 28),
    technical_replicates = list(), expected_final_pairs = 16
  ),
  GSE25070 = list(
    gse = "GSE25070", gpl = "GPL6883", pattern = "block_tumor_normal",
    n_pairs = 26, sample_subset = NULL,
    drop_indices = c(5, 11, 16, 18, 31, 37, 42, 44),
    technical_replicates = list(), expected_final_pairs = 22
  ),
  GSE21510 = list(
    gse = "GSE21510", gpl = "GPL570", pattern = "alt_normal_tumor",
    n_pairs = 23, sample_subset = GSE21510_SAMPLES,
    drop_indices = c(11, 12, 15, 16, 25, 26),
    technical_replicates = list(), expected_final_pairs = 20
  )
)


###############################################################################
# 4. GEO INPUT + PLATFORM ANNOTATION
#
# Local GEO Series Matrix and GPL annotation files are used when available.
# Network download is used only as a fallback.
###############################################################################

read_local_series_matrix <- function(path) {
  message("Using local GEO matrix: ", path)
  
  x <- data.table::fread(
    path,
    skip = "!series_matrix_table_begin",
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE
  )
  
  if (ncol(x) < 2) {
    stop("Could not parse GEO Series Matrix file: ", path)
  }
  
  id_col <- colnames(x)[1]
  x <- x[!grepl("^!series_matrix_table_end", as.character(x[[id_col]])), , drop = FALSE]
  
  probe_id <- gsub('"', "", as.character(x[[id_col]]), fixed = TRUE)
  expr <- x[, -1, drop = FALSE]
  
  # fread normally removes quotes from column names, but strip any remaining ones.
  colnames(expr) <- gsub('"', "", colnames(expr), fixed = TRUE)
  
  expr <- as.data.frame(lapply(expr, function(z) suppressWarnings(as.numeric(z))),
                        check.names = FALSE)
  rownames(expr) <- probe_id
  as.matrix(expr)
}

download_geo_matrix_fallback <- function(gse) {
  # Direct HTTPS fallback for the standard single-platform Series Matrix name.
  # Most of the cohorts used here have a standard matrix file.
  prefix <- sub("[0-9]{3}$", "nnn", gse)
  url <- sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
    prefix, gse, gse
  )
  dest <- file.path(DIR_GEO, paste0(gse, "_series_matrix.txt.gz"))
  
  message("Trying direct HTTPS GEO download: ", url)
  tryCatch(
    {
      utils::download.file(url, destfile = dest, mode = "wb", method = "libcurl", quiet = FALSE)
      if (!file.exists(dest) || file.info(dest)$size <= 0) {
        stop("Downloaded file is empty.")
      }
      dest
    },
    error = function(e) {
      if (file.exists(dest)) unlink(dest)
      stop(
        gse, ": GEO download failed. If you already have the original file, place it at:\n  ",
        file.path(DIR_DATA, paste0(gse, ".gz")),
        "\nOriginal download error: ", conditionMessage(e)
      )
    }
  )
}

load_expression_matrix <- function(cfg) {
  # 1) Best option: use the exact local file used by the original analysis.
  local_candidates <- c(
    file.path(DIR_DATA, paste0(cfg$gse, ".gz")),
    file.path(DIR_DATA, paste0(cfg$gse, "_series_matrix.txt.gz")),
    file.path(DIR_GEO, paste0(cfg$gse, "_series_matrix.txt.gz"))
  )
  local_candidates <- local_candidates[file.exists(local_candidates)]
  
  if (length(local_candidates) > 0) {
    return(list(
      expr = read_local_series_matrix(local_candidates[1]),
      eset = NULL,
      source = local_candidates[1]
    ))
  }
  
  # 2) Try GEOquery without AnnotGPL/getGPL.  AnnotGPL is deliberately FALSE:
  #    annotated GPL downloads are a common source of HTTP 403 errors.
  message("No local matrix found; trying GEOquery without AnnotGPL...")
  gse_obj <- tryCatch(
    GEOquery::getGEO(
      cfg$gse,
      GSEMatrix = TRUE,
      AnnotGPL = FALSE,
      getGPL = FALSE,
      destdir = DIR_GEO
    ),
    error = function(e) e
  )
  
  if (!inherits(gse_obj, "error")) {
    if (!is.list(gse_obj)) gse_obj <- list(gse_obj)
    ann <- vapply(gse_obj, Biobase::annotation, character(1))
    hit <- which(ann == cfg$gpl)
    
    # Some Series Matrix objects have empty annotation labels.  If there is
    # exactly one object, it is safe to use it for these predefined cohorts.
    if (length(hit) == 0 && length(gse_obj) == 1) hit <- 1
    
    if (length(hit) == 1) {
      eset <- gse_obj[[hit]]
      return(list(
        expr = Biobase::exprs(eset),
        eset = eset,
        source = paste0("GEOquery:", cfg$gse)
      ))
    }
  } else {
    message("GEOquery failed: ", conditionMessage(gse_obj))
  }
  
  # 3) Final network fallback: direct HTTPS Series Matrix download.
  direct_file <- download_geo_matrix_fallback(cfg$gse)
  list(
    expr = read_local_series_matrix(direct_file),
    eset = NULL,
    source = direct_file
  )
}

parse_platform_annotation_df <- function(fd, gpl, probe_ids = NULL) {
  # Harmonize platform row identifiers.
  if ("ID" %in% colnames(fd)) {
    rownames(fd) <- as.character(fd$ID)
  } else if ("ID_REF" %in% colnames(fd)) {
    rownames(fd) <- as.character(fd$ID_REF)
  }
  
  if (!is.null(probe_ids)) {
    missing_probe <- setdiff(probe_ids, rownames(fd))
    if (length(missing_probe) > 0) {
      warning(
        gpl, ": ", length(missing_probe),
        " expression probes were absent from the platform table; they will remain unmapped."
      )
    }
    idx <- match(probe_ids, rownames(fd))
    fd <- fd[idx, , drop = FALSE]
    rownames(fd) <- probe_ids
  }
  
  direct_map <- list(
    GPL20115 = list(
      symbol = c("GeneSymbol", "Gene.Symbol", "GENE_SYMBOL"),
      entrez = c("EntrezGeneID", "ENTREZ_GENE_ID", "Entrez_Gene_ID")
    ),
    GPL21185 = list(
      symbol = c("GENE_SYMBOL", "Gene.Symbol"),
      entrez = c("LOCUSLINK_ID", "ENTREZ_GENE_ID")
    ),
    GPL16699 = list(
      symbol = c("GENE_SYMBOL", "Gene.Symbol"),
      entrez = c("LOCUSLINK_ID", "ENTREZ_GENE_ID")
    ),
    GPL14550 = list(
      symbol = c("GENE_SYMBOL", "Gene.Symbol"),
      entrez = c("GENE", "ENTREZ_GENE_ID")
    ),
    GPL6104 = list(
      symbol = c("Symbol", "Gene.Symbol"),
      entrez = c("Entrez_Gene_ID", "ENTREZ_GENE_ID")
    ),
    GPL570 = list(
      symbol = c("Gene.Symbol", "Gene Symbol", "GENE_SYMBOL"),
      entrez = c("ENTREZ_GENE_ID", "Entrez_Gene_ID")
    ),
    GPL6883 = list(
      symbol = c("Symbol", "Gene.Symbol"),
      entrez = c("Entrez_Gene_ID", "ENTREZ_GENE_ID")
    )
  )
  
  if (gpl == "GPL17586") {
    assignment_col <- find_annotation_col(fd, c("gene_assignment", "Gene.Assignment"))
    if (is.na(assignment_col)) {
      stop("GPL17586: gene_assignment column not found.")
    }
    
    # GPL17586 gene_assignment entries are separated by " // "; the fifth field is the first Entrez assignment used by the original pipeline.
    parts <- strsplit(as.character(fd[[assignment_col]]), " // ", fixed = TRUE)
    symbol <- vapply(parts, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
    entrez <- vapply(parts, function(z) if (length(z) >= 5) z[5] else NA_character_, character(1))
  } else {
    mp <- direct_map[[gpl]]
    if (is.null(mp)) stop("No annotation parser configured for ", gpl)
    
    scol <- find_annotation_col(fd, mp$symbol)
    ecol <- find_annotation_col(fd, mp$entrez)
    
    if (is.na(scol) && is.na(ecol)) {
      stop(
        gpl, ": neither symbol nor Entrez annotation was found.\nAvailable columns: ",
        paste(colnames(fd), collapse = ", ")
      )
    }
    
    symbol <- if (!is.na(scol)) fd[[scol]] else rep(NA_character_, nrow(fd))
    entrez <- if (!is.na(ecol)) fd[[ecol]] else rep(NA_character_, nrow(fd))
  }
  
  symbol <- first_symbol(symbol)
  entrez <- first_entrez(entrez)
  
  n_mapped <- sum(!is.na(entrez))
  n_unique <- length(unique(entrez[!is.na(entrez)]))
  
  message(
    gpl, ": mapped ", n_mapped, " / ", length(entrez),
    " expression probes using the platform Entrez/LocusLink annotation; ",
    n_unique, " unique Entrez genes."
  )
  
  if (n_mapped == 0 || n_unique < 1000) {
    stop(
      gpl, ": annotation produced only ", n_unique, " unique Entrez genes.\n",
      "This is too low for a whole-genome expression platform and indicates an ",
      "annotation/probe-ID mismatch.\nAvailable GPL columns: ",
      paste(colnames(fd), collapse = ", ")
    )
  }
  
  tibble(
    probe = if (is.null(probe_ids)) rownames(fd) else probe_ids,
    symbol = symbol,
    entrez = entrez
  )
}

extract_platform_annotation <- function(eset = NULL, gpl, probe_ids) {
  # Prefer the exact GPL text file already used by the original scripts.
  local_gpl <- file.path(DIR_DATA, paste0(gpl, ".txt"))
  
  if (file.exists(local_gpl)) {
    message("Using local GPL annotation: ", local_gpl)
    fd <- utils::read.delim(
      local_gpl,
      comment.char = "#",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      quote = ""
    )
    return(parse_platform_annotation_df(fd, gpl, probe_ids))
  }
  
  # If the Series Matrix ExpressionSet already contains useful feature data,
  # use it without another network request.
  if (!is.null(eset)) {
    fd <- Biobase::fData(eset)
    if (nrow(fd) > 0) {
      try_from_eset <- tryCatch(
        parse_platform_annotation_df(fd, gpl, probe_ids),
        error = function(e) e
      )
      if (!inherits(try_from_eset, "error")) return(try_from_eset)
    }
  }
  
  stop(
    gpl, ": no usable local GPL annotation was found.\n",
    "Place the original platform annotation file here:\n  ",
    local_gpl,
    "\nThis is preferable to redownloading and matches your original analysis."
  )
}

###############################################################################
# 5. LOAD + PREPROCESS ONE GEO COHORT
###############################################################################

load_geo_cohort <- function(cfg) {
  message("\n===== Loading ", cfg$gse, " (", cfg$gpl, ") =====")
  
  loaded <- load_expression_matrix(cfg)
  expr <- loaded$expr
  eset <- loaded$eset
  message("Expression source: ", loaded$source)
  
  # Reproduce the original manually selected sample subset where applicable.
  if (!is.null(cfg$sample_subset)) {
    missing_samples <- setdiff(cfg$sample_subset, colnames(expr))
    if (length(missing_samples) > 0) {
      stop(cfg$gse, ": missing configured samples: ", paste(missing_samples, collapse = ", "))
    }
    expr <- expr[, cfg$sample_subset, drop = FALSE]
  }
  
  # Average the three technical-replicate pairs in GSE84984 exactly as before.
  if (length(cfg$technical_replicates) > 0) {
    for (pair in cfg$technical_replicates) {
      if (!all(pair %in% colnames(expr))) {
        stop(cfg$gse, ": technical replicate(s) not found: ", paste(pair, collapse = ", "))
      }
      expr[, pair[1]] <- rowMeans(expr[, pair, drop = FALSE], na.rm = TRUE)
      expr <- expr[, setdiff(colnames(expr), pair[-1]), drop = FALSE]
    }
  }
  
  # GEO processed matrices are often already log2-scaled. Transform only if needed.
  expr <- log2_if_needed(expr)
  
  expected_before_drop <- 2 * cfg$n_pairs
  if (ncol(expr) != expected_before_drop) {
    stop(
      cfg$gse, ": sample count mismatch after configured subset/replicate handling. ",
      "Expected ", expected_before_drop, ", observed ", ncol(expr), "."
    )
  }
  
  condition <- make_condition(cfg$pattern, cfg$n_pairs)
  patient_local <- make_patient_id(cfg$pattern, cfg$n_pairs)
  
  meta <- tibble(
    sample = colnames(expr),
    condition = condition,
    patient_local = patient_local
  )
  
  # Preserve the exact pair exclusions used in the original pipeline.
  if (length(cfg$drop_indices) > 0) {
    keep <- setdiff(seq_len(ncol(expr)), cfg$drop_indices)
    expr <- expr[, keep, drop = FALSE]
    meta <- meta[keep, , drop = FALSE]
  }
  
  # Renumber retained patients consecutively within cohort.
  old_patients <- unique(meta$patient_local)
  pat_map <- setNames(seq_along(old_patients), old_patients)
  meta$patient_local <- unname(pat_map[as.character(meta$patient_local)])
  meta$patient <- paste0(cfg$gse, "_P", sprintf("%03d", meta$patient_local))
  meta$cohort <- cfg$gse
  meta$condition <- factor(meta$condition, levels = c("normal", "tumor"))
  
  check_pairs(meta, cfg$gse)
  n_pairs_final <- dplyr::n_distinct(meta$patient_local)
  if (n_pairs_final != cfg$expected_final_pairs) {
    stop(cfg$gse, ": expected ", cfg$expected_final_pairs,
         " retained pairs; observed ", n_pairs_final)
  }
  
  # Probe -> Entrez annotation and probe aggregation.
  anno <- extract_platform_annotation(
    eset = eset,
    gpl = cfg$gpl,
    probe_ids = rownames(expr)
  )
  expr_gene <- aggregate_probes_to_entrez(expr, anno$entrez)
  
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expr_gene),
    keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
  )
  
  # Publication-quality cohort PCA (top variable genes only).
  n_top <- min(1000L, nrow(expr_gene))
  rv <- matrixStats::rowVars(expr_gene, na.rm = TRUE)
  top <- order(rv, decreasing = TRUE)[seq_len(n_top)]
  pca <- prcomp(t(expr_gene[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  var_exp <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
  pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("sample") |>
    left_join(meta, by = "sample")
  
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, group = patient)) +
    geom_line(color = "grey75", linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 2.3, alpha = 0.9) +
    labs(
      title = cfg$gse,
      subtitle = paste0(n_pairs_final, " paired patients"),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2]),
      color = "Tissue"
    )
  save_plot(p_pca, paste0("QC_PCA_", cfg$gse), 5.5, 4.5)
  
  # Sample-distribution plot.
  box_df <- as.data.frame(expr_gene) |>
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) |>
    pivot_longer(everything(), names_to = "sample", values_to = "median_expression") |>
    left_join(meta, by = "sample")
  p_box <- ggplot(box_df, aes(x = reorder(sample, median_expression), y = median_expression,
                              fill = condition)) +
    geom_col(width = 0.8) +
    coord_flip() +
    labs(title = paste0(cfg$gse, " sample medians"), x = NULL,
         y = "Median log2 expression", fill = "Tissue") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  save_plot(p_box, paste0("QC_Medians_", cfg$gse), 5.5, 6.5)
  
  list(
    gse = cfg$gse,
    gpl = cfg$gpl,
    expr = expr_gene,
    symbol = symbols,
    meta = meta
  )
}

###############################################################################
# 6. LOAD ALL NINE GEO COHORTS
###############################################################################

geo_data <- purrr::map(geo_cfg, load_geo_cohort)

cohort_summary <- purrr::map_dfr(geo_data, function(x) {
  tibble(
    cohort = x$gse,
    n_pairs = dplyr::n_distinct(x$meta$patient),
    n_samples = nrow(x$meta),
    n_genes = nrow(x$expr)
  )
})

message("Paired GEO patients retained: ", sum(cohort_summary$n_pairs))
write_tsv_safe(cohort_summary, "GEO_cohort_summary.tsv")

p_cohort_summary <- ggplot(
  cohort_summary,
  aes(
    x = n_pairs,
    y = reorder(
      cohort,
      n_pairs
    )
  )
) +
  geom_col(
    fill = unname(
      PLOT_COLORS["down"]
    ),
    width = 0.72
  ) +
  labs(
    title = "Paired patients retained from each GEO cohort",
    x = "Paired patients",
    y = NULL
  )

save_supplementary_plot(
  p_cohort_summary,
  "Figure_S00_GEO_cohort_sizes",
  width = 7,
  height = 5.5
)


###############################################################################
# 7. CORRECTED PAIRED LIMMA WITHIN EACH COHORT
#    Differential expression uses the cohort-scale log2 values.
###############################################################################

run_paired_limma <- function(x) {
  meta <- x$meta
  expr <- x$expr[, meta$sample, drop = FALSE]
  meta$patient <- factor(meta$patient)
  meta$condition <- factor(meta$condition, levels = c("normal", "tumor"))
  
  design <- model.matrix(~ patient + condition, data = meta)
  coef_name <- "conditiontumor"
  if (!coef_name %in% colnames(design)) stop(x$gse, ": condition coefficient missing.")
  
  fit0 <- limma::lmFit(expr, design)
  fit <- limma::eBayes(fit0, trend = TRUE)
  tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  
  # Exact moderated standard error used by limma's moderated t statistic.
  coef_index <- match(coef_name, colnames(design))
  gene_index <- match(rownames(tt), rownames(fit$coefficients))
  se_mod <- fit$stdev.unscaled[gene_index, coef_index] * sqrt(fit$s2.post[gene_index])
  se_mod <- as.numeric(se_mod)
  se_mod[!is.finite(se_mod) | se_mod <= 0] <- NA_real_
  
  tibble(
    entrez = rownames(tt),
    cohort = x$gse,
    n_pairs = dplyr::n_distinct(meta$patient),
    logFC = tt$logFC,
    SE = se_mod,
    AveExpr = tt$AveExpr,
    t = tt$t,
    P.Value = tt$P.Value,
    adj.P.Val = tt$adj.P.Val
  )
}

limma_by_cohort <- purrr::map(geo_data, run_paired_limma)
names(limma_by_cohort) <- names(geo_data)

walk2(limma_by_cohort, names(limma_by_cohort), ~
        write_tsv_safe(.x, paste0("limma_paired_", .y, ".tsv")))

limma_long <- bind_rows(limma_by_cohort)
write_tsv_safe(limma_long, "limma_all_cohorts_long.tsv")

###############################################################################
# 8. RANDOM-EFFECTS META-ANALYSIS ACROSS COHORTS
###############################################################################

safe_rma <- function(df) {
  df <- df |>
    filter(is.finite(logFC), is.finite(SE), SE > 0)
  
  if (nrow(df) < META_MIN_STUDIES) {
    return(tibble(
      k = nrow(df), meta_logFC = NA_real_, meta_SE = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, I2 = NA_real_
    ))
  }
  
  fit <- tryCatch(
    metafor::rma.uni(yi = df$logFC, sei = df$SE, method = META_METHOD),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      k = nrow(df), meta_logFC = NA_real_, meta_SE = NA_real_,
      ci_lb = NA_real_, ci_ub = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, I2 = NA_real_
    ))
  }
  
  tibble(
    k = fit$k,
    meta_logFC = as.numeric(fit$b[1, 1]),
    meta_SE = fit$se,
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    p_value = fit$pval,
    tau2 = fit$tau2,
    I2 = fit$I2
  )
}

message("\nRunning random-effects meta-analysis across genes ...")
meta_results <- limma_long |>
  group_by(entrez) |>
  group_modify(~ safe_rma(.x)) |>
  ungroup()

mean_expr <- limma_long |>
  group_by(entrez) |>
  summarise(
    mean_AveExpr = weighted.mean(AveExpr, w = n_pairs, na.rm = TRUE),
    .groups = "drop"
  )

meta_results <- meta_results |>
  left_join(mean_expr, by = "entrez") |>
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    symbol = unname(AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = entrez, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
    )),
    direction = case_when(
      FDR < META_FDR_CUTOFF & meta_logFC >=  META_LOGFC_CUTOFF ~ "Up",
      FDR < META_FDR_CUTOFF & meta_logFC <= -META_LOGFC_CUTOFF ~ "Down",
      TRUE ~ "Not significant"
    )
  ) |>
  arrange(FDR)

write_tsv_safe(meta_results, "META_all_genes.tsv")

meta_up <- meta_results |> filter(direction == "Up")
meta_down <- meta_results |> filter(direction == "Down")
meta_deg <- meta_results |> filter(direction %in% c("Up", "Down"))

write_tsv_safe(meta_up, "META_DEG_up.tsv")
write_tsv_safe(meta_down, "META_DEG_down.tsv")
write_tsv_safe(meta_deg, "META_DEG_all.tsv")

message("Corrected meta-analysis DEGs: ", nrow(meta_deg),
        " (up=", nrow(meta_up), ", down=", nrow(meta_down), ")")

###############################################################################
# 9. FIGURE 1 — META VOLCANO + MA
###############################################################################

plot_meta <- meta_results |>
  filter(is.finite(meta_logFC), is.finite(FDR), is.finite(mean_AveExpr)) |>
  mutate(
    minus_log10_fdr = -log10(pmax(FDR, .Machine$double.xmin)),
    label = ifelse(symbol %in% FOCAL_GENES, symbol, NA_character_)
  )

p_volcano <- ggplot(plot_meta, aes(meta_logFC, minus_log10_fdr)) +
  geom_point(aes(color = direction), alpha = 0.55, size = 1.2) +
  geom_vline(xintercept = c(-META_LOGFC_CUTOFF, META_LOGFC_CUTOFF),
             linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(META_FDR_CUTOFF),
             linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_label_repel(
    data = subset(plot_meta, !is.na(label)), aes(label = label),
    size = 3, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = c("Down" = "#2C7BB6", "Not significant" = "grey78", "Up" = "#D7191C")) +
  labs(
    title = "A  Meta-analysis volcano plot",
    x = expression("Pooled log"[2]*" fold change (tumor vs normal)"),
    y = expression(-log[10](FDR)), color = NULL
  )

p_ma <- ggplot(plot_meta, aes(mean_AveExpr, meta_logFC)) +
  geom_point(aes(color = direction), alpha = 0.55, size = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-META_LOGFC_CUTOFF, META_LOGFC_CUTOFF),
             linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_label_repel(
    data = subset(plot_meta, !is.na(label)), aes(label = label),
    size = 3, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = c("Down" = "#2C7BB6", "Not significant" = "grey78", "Up" = "#D7191C")) +
  labs(
    title = "B  Meta-analysis MA plot",
    x = "Mean log2 expression across cohorts",
    y = expression("Pooled log"[2]*" fold change"), color = NULL
  ) +
  theme(legend.position = "none")

fig1 <- p_volcano + p_ma + patchwork::plot_layout(guides = "collect")

###############################################################################
# 10. FOCAL-GENE META-ANALYSIS AND SENSITIVITY PLOTS
###############################################################################

meta_focal_summary <- list()
meta_forest_plots <- list()
leave_one_out_plots <- list()

for (gene in FOCAL_GENES) {
  
  entrez_gene <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = gene,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  ) |>
    unname()
  
  df <- limma_long |>
    filter(
      entrez == entrez_gene,
      is.finite(logFC),
      is.finite(SE),
      SE > 0
    )
  
  if (nrow(df) < META_MIN_STUDIES) {
    next
  }
  
  fit <- metafor::rma.uni(
    yi = df$logFC,
    sei = df$SE,
    method = META_METHOD,
    slab = df$cohort
  )
  
  study_forest <- df |>
    transmute(
      cohort = cohort,
      estimate = logFC,
      ci_low = logFC - 1.96 * SE,
      ci_high = logFC + 1.96 * SE,
      type = "Cohort"
    )
  
  pooled_forest <- tibble(
    cohort = "Random-effects pooled",
    estimate = as.numeric(fit$b[1, 1]),
    ci_low = fit$ci.lb,
    ci_high = fit$ci.ub,
    type = "Pooled"
  )
  
  forest_df <- bind_rows(
    study_forest,
    pooled_forest
  ) |>
    mutate(
      cohort = factor(
        cohort,
        levels = rev(c(as.character(study_forest$cohort), "Random-effects pooled"))
      )
    )
  
  p_forest <- ggplot(
    forest_df,
    aes(
      x = estimate,
      y = cohort,
      color = type
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.45,
      color = "grey55"
    ) +
    geom_errorbarh(
      aes(
        xmin = ci_low,
        xmax = ci_high
      ),
      height = 0.16,
      linewidth = 0.7
    ) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c(
        Cohort = unname(PLOT_COLORS["down"]),
        Pooled = unname(PLOT_COLORS["accent"])
      )
    ) +
    labs(
      title = gene,
      x = expression("log"[2]*" fold change (tumor vs normal)"),
      y = NULL,
      color = NULL
    ) +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 9)
    )
  
  meta_forest_plots[[gene]] <- p_forest
  save_plot(
    p_forest,
    paste0("Forest_", gene),
    width = 7,
    height = 5.5
  )
  
  loo <- metafor::leave1out(fit) |>
    as.data.frame() |>
    tibble::rownames_to_column("omitted_study")
  
  write_tsv_safe(
    loo,
    paste0("Leave_one_out_", gene, ".tsv")
  )
  
  if (
    all(
      c("estimate", "ci.lb", "ci.ub") %in%
      colnames(loo)
    )
  ) {
    
    loo_plot_df <- loo |>
      mutate(
        omitted_study = factor(
          omitted_study,
          levels = rev(omitted_study)
        )
      )
    
    p_loo <- ggplot(
      loo_plot_df,
      aes(
        x = estimate,
        y = omitted_study
      )
    ) +
      geom_vline(
        xintercept = as.numeric(fit$b[1, 1]),
        linetype = "dotted",
        linewidth = 0.6,
        color = PLOT_COLORS["accent"]
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.45,
        color = "grey55"
      ) +
      geom_errorbarh(
        aes(
          xmin = ci.lb,
          xmax = ci.ub
        ),
        height = 0.16,
        linewidth = 0.7,
        color = PLOT_COLORS["down"]
      ) +
      geom_point(
        size = 2.6,
        color = PLOT_COLORS["down"]
      ) +
      labs(
        title = gene,
        subtitle = "Leave-one-study-out random-effects estimate",
        x = expression("Pooled log"[2]*" fold change"),
        y = "Omitted cohort"
      ) +
      theme(
        axis.text.y = element_text(size = 9)
      )
    
    leave_one_out_plots[[gene]] <- p_loo
  }
  
  meta_focal_summary[[gene]] <- tibble(
    gene = gene,
    k = fit$k,
    pooled_logFC = as.numeric(fit$b[1, 1]),
    ci_lb = fit$ci.lb,
    ci_ub = fit$ci.ub,
    p_value = fit$pval,
    I2 = fit$I2,
    tau2 = fit$tau2
  )
}

meta_focal_summary <- bind_rows(meta_focal_summary)

write_tsv_safe(
  meta_focal_summary,
  "Focal_genes_meta_summary.tsv"
)

if (length(meta_forest_plots) > 0) {
  p_meta_forest_supp <- patchwork::wrap_plots(
    meta_forest_plots[intersect(FOCAL_GENES, names(meta_forest_plots))],
    ncol = 3
  ) +
    patchwork::plot_annotation(tag_levels = "A")
  
  save_supplementary_plot(
    p_meta_forest_supp,
    "Figure_S03_Focal_gene_meta_forest",
    width = 14,
    height = 5.5
  )
}

if (length(leave_one_out_plots) > 0) {
  p_loo_supp <- patchwork::wrap_plots(
    leave_one_out_plots[intersect(FOCAL_GENES, names(leave_one_out_plots))],
    ncol = 3
  ) +
    patchwork::plot_annotation(tag_levels = "A")
  
  save_supplementary_plot(
    p_loo_supp,
    "Figure_S04_Leave_one_out_sensitivity",
    width = 14,
    height = 5.5
  )
}

###############################################################################
# 11. CROSS-PLATFORM INTEGRATION FOR WGCNA ONLY
#     Gene-wise Z-score is used here, not for fold-change estimation.
###############################################################################

# ---- Cross-platform gene-overlap QC -----------------------------------------
# The original manuscript reported 15,890 intersecting genes across the final
# nine cohorts.  We do not force that exact historical number (annotations can
# change), but an intersection near zero is biologically impossible and signals
# an annotation-ID problem.  Diagnose it *before* ComBat/WGCNA.
entrez_sets <- lapply(geo_data, function(x) unique(rownames(x$expr)))
names(entrez_sets) <- names(geo_data)

# Per-cohort valid gene counts.
overlap_counts <- tibble(
  cohort = names(entrez_sets),
  n_valid_entrez_genes = vapply(entrez_sets, length, integer(1))
)
write_tsv_safe(overlap_counts, "WGCNA_gene_counts_by_cohort.tsv")

# Progressive intersection in the configured cohort order.
progressive <- vector("list", length(entrez_sets))
running <- NULL
for (i in seq_along(entrez_sets)) {
  running <- if (is.null(running)) entrez_sets[[i]] else intersect(running, entrez_sets[[i]])
  progressive[[i]] <- tibble(
    step = i,
    added_cohort = names(entrez_sets)[i],
    n_common_genes = length(running)
  )
  message("WGCNA overlap after ", names(entrez_sets)[i], ": ", length(running), " genes")
}
progressive_df <- bind_rows(progressive)
write_tsv_safe(progressive_df, "WGCNA_progressive_gene_overlap.tsv")

# Pairwise overlap matrix is useful if one platform is still problematic.
pairwise_overlap <- outer(
  seq_along(entrez_sets), seq_along(entrez_sets),
  Vectorize(function(i, j) length(intersect(entrez_sets[[i]], entrez_sets[[j]])))
)
dimnames(pairwise_overlap) <- list(names(entrez_sets), names(entrez_sets))
write_matrix_tsv(pairwise_overlap, "WGCNA_pairwise_gene_overlap.tsv")

p_overlap_progressive <- ggplot(
  progressive_df,
  aes(
    x = step,
    y = n_common_genes
  )
) +
  geom_line(
    linewidth = 0.8,
    color = PLOT_COLORS["accent"]
  ) +
  geom_point(
    size = 2.5,
    color = PLOT_COLORS["accent"]
  ) +
  scale_x_continuous(
    breaks = progressive_df$step,
    labels = progressive_df$added_cohort
  ) +
  labs(
    title = "Progressive cross-platform gene intersection",
    x = NULL,
    y = "Common genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

pairwise_overlap_long <- as.data.frame(
  as.table(pairwise_overlap),
  stringsAsFactors = FALSE
)

if (ncol(pairwise_overlap_long) != 3L) {
  stop(
    "Unexpected pairwise-overlap table structure: expected 3 columns, observed ",
    ncol(pairwise_overlap_long),
    "."
  )
}

colnames(pairwise_overlap_long) <- c(
  "cohort_1",
  "cohort_2",
  "n_common_genes"
)

pairwise_overlap_long <- tibble::as_tibble(
  pairwise_overlap_long
)

p_overlap_pairwise <- ggplot(
  pairwise_overlap_long,
  aes(
    x = cohort_1,
    y = cohort_2,
    fill = n_common_genes
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = n_common_genes),
    size = 2.7
  ) +
  scale_fill_gradient(
    low = "white",
    high = PLOT_COLORS["down"]
  ) +
  labs(
    title = "Pairwise platform gene overlap",
    x = NULL,
    y = NULL,
    fill = "Genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_supplementary_plot(
  p_overlap_progressive +
    p_overlap_pairwise +
    patchwork::plot_annotation(tag_levels = "A"),
  "Figure_S01_WGCNA_gene_overlap_QC",
  width = 13,
  height = 6
)


common_entrez <- Reduce(intersect, entrez_sets)
message("Common platform-annotated Entrez genes across all 9 WGCNA cohorts: ", length(common_entrez))

# The manuscript dataset contained 15,890 genes in the nine-cohort intersection.
EXPECTED_COMMON_GENES <- 15890L
if (length(common_entrez) != EXPECTED_COMMON_GENES) {
  warning(
    "The current platform intersection contains ",
    length(common_entrez),
    " genes; the manuscript dataset contained 15,890. ",
    "Check the platform annotations before comparing exact WGCNA results."
  )
} else {
  message("Cross-platform intersection: 15,890 genes.")
}

# A catastrophic overlap is almost always annotation corruption.  Stop with a
# diagnostic table rather than allowing ComBat to fail cryptically.
MIN_WGCNA_COMMON_GENES <- 5000L
if (length(common_entrez) < MIN_WGCNA_COMMON_GENES) {
  stop(
    "Only ", length(common_entrez), " platform-annotated Entrez genes are common to all 9 cohorts. ",
    "This is too low for the intended cross-platform WGCNA and indicates a remaining ",
    "platform-annotation mismatch. Check Results_manuscript/Tables/",
    "WGCNA_progressive_gene_overlap.tsv and WGCNA_pairwise_gene_overlap.tsv to identify ",
    "the offending cohort. (The historical pipeline reported 15,890 intersecting genes.)"
  )
}

wgcna_parts <- list()
wgcna_meta <- list()

for (i in seq_along(geo_data)) {
  x <- geo_data[[i]]
  mat <- x$expr[common_entrez, x$meta$sample, drop = FALSE]
  
  # Correct standardization: each gene is standardized across samples within cohort.
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  
  wgcna_parts[[i]] <- z
  wgcna_meta[[i]] <- x$meta
}

merged_z <- do.call(cbind, wgcna_parts)
meta_all <- bind_rows(wgcna_meta)
stopifnot(identical(colnames(merged_z), meta_all$sample))

# Give retained paired patients a global 1..145 index in first-appearance order.
patient_order <- unique(meta_all$patient)
patient_global_map <- setNames(seq_along(patient_order), patient_order)
meta_all$global_patient <- unname(patient_global_map[meta_all$patient])

batch <- factor(meta_all$cohort)
combat_mod <- model.matrix(~ condition, data = meta_all)
combat_z <- sva::ComBat(
  dat = merged_z,
  batch = batch,
  mod = combat_mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

write_matrix_tsv(combat_z, "WGCNA_integrated_ComBat_geneZ.tsv")
write_tsv_safe(meta_all, "WGCNA_sample_metadata_before_WGCNA_outlier_removal.tsv")

# PCA before/after ComBat
pca_from_matrix <- function(mat, meta, label) {
  rv <- matrixStats::rowVars(mat)
  top <- order(rv, decreasing = TRUE)[seq_len(min(1500, length(rv)))]
  pc <- prcomp(t(mat[top, , drop = FALSE]), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  df <- as.data.frame(pc$x[, 1:2, drop = FALSE]) |>
    rownames_to_column("sample") |>
    left_join(meta, by = "sample")
  ggplot(df, aes(PC1, PC2, color = cohort, shape = condition)) +
    geom_point(size = 2, alpha = 0.85) +
    labs(title = label, x = sprintf("PC1 (%.1f%%)", ve[1]),
         y = sprintf("PC2 (%.1f%%)", ve[2]), color = "Cohort", shape = "Tissue")
}

p_before <- pca_from_matrix(merged_z, meta_all, "A  Before ComBat")
p_after  <- pca_from_matrix(combat_z, meta_all, "B  After ComBat")
fig_batch <- p_before + p_after + patchwork::plot_layout(guides = "collect")
save_plot(fig_batch, "QC_Batch_correction_PCA", 12, 5.5)
save_supplementary_plot(
  fig_batch,
  "Figure_S13_Batch_correction_PCA",
  width = 12,
  height = 5.5
)


###############################################################################
# 12. DEG HEATMAP
###############################################################################

heat_entrez <- intersect(
  meta_deg$entrez,
  rownames(combat_z)
)

if (length(heat_entrez) > 1) {
  
  ann_col <- data.frame(
    Tissue = meta_all$condition,
    Cohort = meta_all$cohort,
    row.names = meta_all$sample
  )
  
  heat_mat <- combat_z[
    heat_entrez,
    ,
    drop = FALSE
  ]
  
  heat_mat <- t(
    scale(
      t(heat_mat)
    )
  )
  
  heat_mat[
    !is.finite(heat_mat)
  ] <- 0
  
  heat_colors <- grDevices::colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(101)
  
  pheatmap::pheatmap(
    heat_mat,
    annotation_col = ann_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = heat_colors,
    main = "Meta-analysis DEGs across the integrated paired CRC cohorts",
    filename = file.path(
      DIR_FIG,
      "DEG_heatmap.pdf"
    ),
    width = 11,
    height = 10
  )
  
  pheatmap::pheatmap(
    heat_mat,
    annotation_col = ann_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = heat_colors,
    main = "Meta-analysis DEGs across the integrated paired CRC cohorts",
    filename = file.path(
      DIR_SUPPLEMENTARY,
      "Figure_S14_DEG_heatmap.pdf"
    ),
    width = 11,
    height = 10
  )
  
  pheatmap::pheatmap(
    heat_mat,
    annotation_col = ann_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = heat_colors,
    main = "Meta-analysis DEGs across the integrated paired CRC cohorts",
    filename = file.path(
      DIR_SUPPLEMENTARY,
      "Figure_S14_DEG_heatmap.png"
    ),
    width = 11,
    height = 10
  )
}

###############################################################################
# 13. ENRICHMENT FUNCTIONS
###############################################################################

get_msig_collection <- function(collection) {
  fn_args <- names(formals(msigdbr::msigdbr))
  if ("collection" %in% fn_args) {
    msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
  } else {
    msigdbr::msigdbr(species = "Homo sapiens", category = collection)
  }
}

msig_term2gene <- function(collection) {
  x <- get_msig_collection(collection)
  gene_col <- if ("ncbi_gene" %in% colnames(x)) "ncbi_gene" else "entrez_gene"
  x |>
    dplyr::select(gs_name, all_of(gene_col)) |>
    dplyr::rename(gene = all_of(gene_col)) |>
    mutate(gene = as.character(gene)) |>
    distinct()
}

run_enrichment_set <- function(entrez_genes, universe_entrez, prefix) {
  entrez_genes <- unique(as.character(entrez_genes))
  universe_entrez <- unique(as.character(universe_entrez))
  
  if (length(entrez_genes) < 5) return(NULL)
  
  go_list <- lapply(c("BP", "MF", "CC"), function(ont) {
    clusterProfiler::enrichGO(
      gene = entrez_genes,
      universe = universe_entrez,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      pAdjustMethod = "BH",
      readable = FALSE
    )
  })
  names(go_list) <- c("BP", "MF", "CC")
  
  kegg <- clusterProfiler::enrichKEGG(
    gene = entrez_genes,
    universe = universe_entrez,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  
  msig_H  <- clusterProfiler::enricher(
    gene = entrez_genes, universe = universe_entrez,
    TERM2GENE = msig_term2gene("H"), pvalueCutoff = 1, qvalueCutoff = 1
  )
  msig_C1 <- clusterProfiler::enricher(
    gene = entrez_genes, universe = universe_entrez,
    TERM2GENE = msig_term2gene("C1"), pvalueCutoff = 1, qvalueCutoff = 1
  )
  msig_C3 <- clusterProfiler::enricher(
    gene = entrez_genes, universe = universe_entrez,
    TERM2GENE = msig_term2gene("C3"), pvalueCutoff = 1, qvalueCutoff = 1
  )
  
  for (ont in names(go_list)) {
    write_tsv_safe(as.data.frame(go_list[[ont]]), paste0(prefix, "_GO_", ont, ".tsv"))
  }
  write_tsv_safe(as.data.frame(kegg), paste0(prefix, "_KEGG.tsv"))
  write_tsv_safe(as.data.frame(msig_H), paste0(prefix, "_Hallmark.tsv"))
  write_tsv_safe(as.data.frame(msig_C1), paste0(prefix, "_Cytoband_C1.tsv"))
  write_tsv_safe(as.data.frame(msig_C3), paste0(prefix, "_Regulatory_C3.tsv"))
  
  list(GO = go_list, KEGG = kegg, H = msig_H, C1 = msig_C1, C3 = msig_C3)
}

enrichment_table_for_plot <- function(
    x,
    collection,
    gene_set,
    top_n = 8
) {
  
  if (is.null(x)) {
    return(tibble())
  }
  
  d <- as.data.frame(x)
  
  if (
    nrow(d) == 0 ||
    !"p.adjust" %in% colnames(d)
  ) {
    return(tibble())
  }
  
  description_col <- if (
    "Description" %in% colnames(d)
  ) {
    "Description"
  } else {
    "ID"
  }
  
  d |>
    filter(
      !is.na(p.adjust),
      is.finite(p.adjust),
      p.adjust < 0.05
    ) |>
    arrange(p.adjust) |>
    slice_head(n = top_n) |>
    transmute(
      gene_set = gene_set,
      collection = collection,
      term = .data[[description_col]],
      gene_ratio = if (
        "GeneRatio" %in% colnames(d)
      ) {
        ratio_to_numeric(GeneRatio)
      } else {
        NA_real_
      },
      count = if (
        "Count" %in% colnames(d)
      ) {
        Count
      } else {
        NA_real_
      },
      FDR = p.adjust
    )
}

collect_enrichment_for_plot <- function(
    enrichment_set,
    gene_set,
    top_n = 8
) {
  
  if (is.null(enrichment_set)) {
    return(tibble())
  }
  
  go_rows <- bind_rows(
    lapply(
      names(enrichment_set$GO),
      function(ont) {
        enrichment_table_for_plot(
          enrichment_set$GO[[ont]],
          paste0("GO ", ont),
          gene_set,
          top_n
        )
      }
    )
  )
  
  bind_rows(
    go_rows,
    enrichment_table_for_plot(
      enrichment_set$KEGG,
      "KEGG",
      gene_set,
      top_n
    ),
    enrichment_table_for_plot(
      enrichment_set$H,
      "Hallmark",
      gene_set,
      top_n
    ),
    enrichment_table_for_plot(
      enrichment_set$C1,
      "C1 cytoband",
      gene_set,
      top_n
    ),
    enrichment_table_for_plot(
      enrichment_set$C3,
      "C3 regulatory",
      gene_set,
      top_n
    )
  )
}

plot_enrichment_collection <- function(
    df,
    title
) {
  
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  ggplot(
    df,
    aes(
      x = gene_ratio,
      y = reorder(term, gene_ratio),
      size = count,
      color = -log10(FDR)
    )
  ) +
    geom_point(alpha = 0.9) +
    facet_grid(
      collection ~ gene_set,
      scales = "free_y",
      space = "free_y"
    ) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      size = "Genes",
      color = expression(-log[10](FDR))
    ) +
    theme(
      axis.text.y = element_text(size = 7)
    )
}


universe_entrez <- meta_results |>
  filter(is.finite(p_value)) |>
  pull(entrez) |>
  unique()

E_UP   <- run_enrichment_set(meta_up$entrez, universe_entrez, "DEG_UP")
E_DOWN <- run_enrichment_set(meta_down$entrez, universe_entrez, "DEG_DOWN")

deg_enrichment_supp <- bind_rows(
  collect_enrichment_for_plot(
    E_UP,
    "Upregulated DEGs",
    top_n = 6
  ),
  collect_enrichment_for_plot(
    E_DOWN,
    "Downregulated DEGs",
    top_n = 6
  )
)

p_deg_enrichment_supp <- plot_enrichment_collection(
  deg_enrichment_supp,
  "Functional enrichment of meta-analysis DEGs"
)

if (!is.null(p_deg_enrichment_supp)) {
  save_supplementary_plot(
    p_deg_enrichment_supp,
    "Figure_S07_DEG_full_enrichment",
    width = 14,
    height = 18
  )
}


###############################################################################
# 14. GO ENRICHMENT PLOT DATA
###############################################################################

extract_go_for_plot <- function(
    enrich_obj,
    direction,
    top_n = 10
) {
  
  dplyr::bind_rows(
    lapply(
      names(
        enrich_obj$GO
      ),
      function(ont) {
        
        as.data.frame(
          enrich_obj$GO[[ont]]
        ) |>
          dplyr::filter(
            p.adjust < 0.05
          ) |>
          dplyr::arrange(
            p.adjust
          ) |>
          dplyr::slice_head(
            n = top_n
          ) |>
          dplyr::mutate(
            ontology = ont,
            direction = direction,
            GeneRatioNum = ratio_to_numeric(
              GeneRatio
            )
          )
      }
    )
  )
}

go_plot_df <- dplyr::bind_rows(
  extract_go_for_plot(
    E_UP,
    "Up"
  ),
  extract_go_for_plot(
    E_DOWN,
    "Down"
  )
)


###############################################################################
# 15. KEGG ENRICHMENT PLOT DATA
###############################################################################

extract_kegg_for_plot <- function(
    enrich_obj,
    direction,
    top_n = 10
) {
  
  as.data.frame(
    enrich_obj$KEGG
  ) |>
    dplyr::filter(
      p.adjust < 0.05
    ) |>
    dplyr::arrange(
      p.adjust
    ) |>
    dplyr::slice_head(
      n = top_n
    ) |>
    dplyr::mutate(
      direction = direction,
      GeneRatioNum = ratio_to_numeric(
        GeneRatio
      )
    )
}

kegg_plot_df <- dplyr::bind_rows(
  extract_kegg_for_plot(
    E_UP,
    "Up"
  ),
  extract_kegg_for_plot(
    E_DOWN,
    "Down"
  )
)


###############################################################################
# 16. WGCNA
###############################################################################

keep_wgcna_samples <- !meta_all$global_patient %in% WGCNA_OUTLIER_GLOBAL_PATIENTS
meta_wgcna <- meta_all[keep_wgcna_samples, , drop = FALSE]
datExpr <- t(combat_z[, meta_wgcna$sample, drop = FALSE])

message("WGCNA samples: ", nrow(datExpr), "; paired patients: ", dplyr::n_distinct(meta_wgcna$patient))

# Use WGCNA correlation functions explicitly during network construction.
cor   <- WGCNA::cor
bicor <- WGCNA::bicor

# Basic WGCNA QC
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  meta_wgcna <- meta_wgcna[rownames(datExpr), , drop = FALSE]
}

# Sample dendrogram
sample_tree <- hclust(
  dist(datExpr),
  method = "average"
)

draw_sample_dendrogram <- function() {
  plot(
    sample_tree,
    main = "WGCNA sample clustering after paired outlier removal",
    xlab = "",
    sub = "",
    cex = 0.5
  )
}

pdf(
  file.path(
    DIR_FIG,
    "WGCNA_sample_dendrogram.pdf"
  ),
  width = 12,
  height = 6
)
draw_sample_dendrogram()
dev.off()

png(
  file.path(
    DIR_SUPPLEMENTARY,
    "Figure_S05_WGCNA_sample_dendrogram.png"
  ),
  width = 4800,
  height = 2400,
  res = 400
)
draw_sample_dendrogram()
dev.off()

pdf(
  file.path(
    DIR_SUPPLEMENTARY,
    "Figure_S05_WGCNA_sample_dendrogram.pdf"
  ),
  width = 12,
  height = 6
)
draw_sample_dendrogram()
dev.off()

powers <- c(1:11, seq(12, 40, 2))
sft <- WGCNA::pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
)
sft_df <- sft$fitIndices

# Select the lowest power that reaches the pre-specified scale-free fit target
# with a negative slope. This avoids forcing a preset power onto the
# corrected integrated matrix if its topology does not support that choice.
eligible_power <- sft_df |>
  dplyr::filter(is.finite(SFT.R.sq), is.finite(slope),
                SFT.R.sq >= WGCNA_SFT_R2_TARGET, slope < 0) |>
  dplyr::arrange(Power)

if (nrow(eligible_power) > 0) {
  WGCNA_SOFT_POWER_USED <- eligible_power$Power[1]
} else {
  negative_slope <- sft_df |>
    dplyr::filter(is.finite(SFT.R.sq), is.finite(slope), slope < 0) |>
    dplyr::arrange(dplyr::desc(SFT.R.sq), Power)
  if (nrow(negative_slope) == 0) {
    stop("No usable soft-threshold power with a negative scale-free fit slope was found.")
  }
  WGCNA_SOFT_POWER_USED <- negative_slope$Power[1]
  warning(
    "No tested power reached R^2 >= ", WGCNA_SFT_R2_TARGET,
    "; using the best negative-slope candidate: power=", WGCNA_SOFT_POWER_USED,
    " (R^2=", signif(negative_slope$SFT.R.sq[1], 3), ")."
  )
}

message(
  "WGCNA soft-threshold selected from corrected data: power=",
  WGCNA_SOFT_POWER_USED
)

sft_df <- sft_df |>
  dplyr::mutate(selected = Power == WGCNA_SOFT_POWER_USED)
write_tsv_safe(sft_df, "WGCNA_soft_threshold_diagnostics.tsv")

# Publication-quality soft-threshold diagnostics are rendered in Section 28.

net <- WGCNA::blockwiseModules(
  datExpr,
  power = WGCNA_SOFT_POWER_USED,
  maxBlockSize = 16000,
  minModuleSize = WGCNA_MIN_MODULE_SIZE,
  mergeCutHeight = WGCNA_MERGE_CUT,
  networkType = "signed",
  TOMType = "signed",
  useCorOptionsThroughout = TRUE,
  numericLabels = FALSE,
  randomSeed = 1234,
  verbose = 5
)

module_colors <- net$colors
MEs <- WGCNA::orderMEs(net$MEs)

# Dendrogram and module colors
draw_wgcna_dendrogram <- function() {
  if (length(net$dendrograms) == 1) {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      cbind(net$unmergedColors, net$colors),
      c("Unmerged", "Merged"),
      dendroLabels = FALSE,
      addGuide = TRUE,
      hang = 0.03,
      guideHang = 0.05
    )
  } else {
    for (b in seq_along(net$dendrograms)) {
      WGCNA::plotDendroAndColors(
        net$dendrograms[[b]],
        net$colors[net$blockGenes[[b]]],
        paste0("Block ", b),
        dendroLabels = FALSE,
        addGuide = TRUE
      )
    }
  }
}

pdf(
  file.path(
    DIR_FIG,
    "Figure_5C_WGCNA_dendrogram.pdf"
  ),
  width = 12,
  height = 6
)
draw_wgcna_dendrogram()
dev.off()

pdf(
  file.path(
    DIR_SUPPLEMENTARY,
    "Figure_S02_WGCNA_dendrogram.pdf"
  ),
  width = 12,
  height = 6
)
draw_wgcna_dendrogram()
dev.off()

png(
  file.path(
    DIR_SUPPLEMENTARY,
    "Figure_S02_WGCNA_dendrogram.png"
  ),
  width = 4800,
  height = 2400,
  res = 400
)
draw_wgcna_dendrogram()
dev.off()

module_map <- tibble(
  entrez = colnames(datExpr),
  module = module_colors,
  symbol = unname(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = colnames(datExpr), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
  ))
)

write_tsv_safe(module_map, "WGCNA_module_membership_map.tsv")
write_tsv_safe(as.data.frame(table(module_colors)), "WGCNA_module_sizes.tsv")

module_size_df <- tibble(
  module = names(table(module_colors)),
  n_genes = as.integer(table(module_colors))
) |>
  arrange(n_genes)

p_module_sizes <- ggplot(
  module_size_df,
  aes(
    x = n_genes,
    y = reorder(module, n_genes),
    fill = module
  )
) +
  geom_col(width = 0.75) +
  scale_fill_identity() +
  labs(
    title = "WGCNA module sizes",
    x = "Genes",
    y = NULL
  )

save_supplementary_plot(
  p_module_sizes,
  "Figure_S06_WGCNA_module_sizes",
  width = 7,
  height = 6
)


###############################################################################
# 17. PAIRED MODULE–TISSUE ASSOCIATION + PAIRED GENE SIGNIFICANCE
###############################################################################

meta_wgcna$patient <- factor(meta_wgcna$patient)
meta_wgcna$condition <- factor(meta_wgcna$condition, levels = c("normal", "tumor"))
design_w <- model.matrix(~ patient + condition, data = meta_wgcna)
coef_w <- which(colnames(design_w) == "conditiontumor")

# Module eigengene association, accounting for patient pairing.
fit_me0 <- limma::lmFit(t(MEs), design_w)
fit_me  <- limma::eBayes(fit_me0, trend = TRUE)
me_tt <- limma::topTable(fit_me, coef = coef_w, number = Inf, sort.by = "none")

me_se_raw <- fit_me0$stdev.unscaled[, coef_w] * fit_me0$sigma
me_t_raw  <- fit_me0$coefficients[, coef_w] / me_se_raw
me_df      <- fit_me0$df.residual
me_partial_r <- me_t_raw / sqrt(me_t_raw^2 + me_df)

module_assoc <- tibble(
  ME = rownames(me_tt),
  module = sub("^ME", "", rownames(me_tt)),
  paired_effect = me_tt$logFC,
  partial_r = unname(me_partial_r[rownames(me_tt)]),
  p_value = me_tt$P.Value,
  FDR = p.adjust(me_tt$P.Value, "BH")
) |>
  arrange(desc(partial_r))

write_tsv_safe(module_assoc, "WGCNA_paired_module_trait_association.tsv")

# Publication-quality module–tissue visualization is rendered in Section 28.

# Gene significance, also accounting for patient pairing.
fit_g0 <- limma::lmFit(t(datExpr), design_w)
fit_g  <- limma::eBayes(fit_g0, trend = TRUE)

g_se_raw <- fit_g0$stdev.unscaled[, coef_w] * fit_g0$sigma
g_t_raw  <- fit_g0$coefficients[, coef_w] / g_se_raw
g_df      <- fit_g0$df.residual
g_partial_r <- g_t_raw / sqrt(g_t_raw^2 + g_df)
names(g_partial_r) <- rownames(fit_g0$coefficients)

# Module membership (standard WGCNA correlation with own module eigengene).
MM <- WGCNA::cor(datExpr, MEs, use = "p")

module_map <- module_map |>
  mutate(
    me_name = paste0("ME", module),
    MM = mapply(function(g, me) {
      if (me %in% colnames(MM) && g %in% rownames(MM)) MM[g, me] else NA_real_
    }, entrez, me_name),
    GS_signed = g_partial_r[entrez],
    GS = abs(GS_signed),
    hub = module != "grey" & abs(MM) >= WGCNA_MM_CUTOFF & GS >= WGCNA_GS_CUTOFF
  ) |>
  dplyr::select(-me_name)

write_tsv_safe(module_map, "WGCNA_gene_MM_GS_hubs.tsv")

# Preserve manuscript color names if they exist; otherwise use strongest modules.
normal_module <- module_assoc$module[which.min(module_assoc$partial_r)]
tumor_module  <- module_assoc$module[which.max(module_assoc$partial_r)]
blue_module <- if ("blue" %in% module_map$module) "blue" else normal_module
turquoise_module <- if ("turquoise" %in% module_map$module) "turquoise" else tumor_module

blue_hubs <- module_map |> filter(module == blue_module, hub)
turquoise_hubs <- module_map |> filter(module == turquoise_module, hub)

write_tsv_safe(blue_hubs, paste0("WGCNA_hubs_", blue_module, ".tsv"))
write_tsv_safe(turquoise_hubs, paste0("WGCNA_hubs_", turquoise_module, ".tsv"))

###############################################################################
# 18. WGCNA HUB AND DEG OVERLAP DATA
###############################################################################

# The final hub-gene scatter plots and DEG-overlap diagrams are rendered in
# Section 28 from module_map, the module-specific hub tables, and meta-analysis
# DEG sets.


###############################################################################
# 19. WGCNA MODULE ENRICHMENT
###############################################################################

E_TURQ <- run_enrichment_set(
  module_map |> filter(module == turquoise_module) |> pull(entrez),
  universe_entrez,
  paste0("WGCNA_", turquoise_module)
)
E_BLUE <- run_enrichment_set(
  module_map |> filter(module == blue_module) |> pull(entrez),
  universe_entrez,
  paste0("WGCNA_", blue_module)
)

wgcna_enrichment_supp <- bind_rows(
  collect_enrichment_for_plot(
    E_TURQ,
    paste0(turquoise_module, " module"),
    top_n = 6
  ),
  collect_enrichment_for_plot(
    E_BLUE,
    paste0(blue_module, " module"),
    top_n = 6
  )
)

p_wgcna_enrichment_supp <- plot_enrichment_collection(
  wgcna_enrichment_supp,
  "Functional enrichment of WGCNA modules"
)

if (!is.null(p_wgcna_enrichment_supp)) {
  save_supplementary_plot(
    p_wgcna_enrichment_supp,
    "Figure_S08_WGCNA_module_enrichment",
    width = 14,
    height = 18
  )
}


###############################################################################
# 20. STRING PPI + CENTRALITY
#     If an exact Cytoscape centrality export is available, place it at:
#       Data/PPI_centrality_export.tsv
#     with a 'symbol' column and the original centrality columns.
###############################################################################

build_string_ppi <- function(deg_symbols) {
  string_db <- STRINGdb::STRINGdb$new(
    version = STRING_VERSION,
    species = 9606,
    score_threshold = STRING_SCORE_THRESHOLD,
    input_directory = DIR_STRING
  )
  
  mapped <- string_db$map(
    data.frame(symbol = unique(na.omit(deg_symbols))),
    "symbol",
    removeUnmappedRows = TRUE
  )
  
  ints <- string_db$get_interactions(mapped$STRING_id)
  id2sym <- setNames(mapped$symbol, mapped$STRING_id)
  
  edges <- tibble(
    from = unname(id2sym[ints$from]),
    to = unname(id2sym[ints$to]),
    combined_score = ints$combined_score
  ) |>
    filter(!is.na(from), !is.na(to), from != to) |>
    distinct(pmin(from, to), pmax(from, to), .keep_all = TRUE) |>
    dplyr::select(from, to, combined_score)
  
  vertices <- tibble(name = sort(unique(c(edges$from, edges$to))))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
  
  list(graph = g, edges = edges, mapping = mapped)
}

ppi <- retry_remote_call(
  label = "STRING PPI retrieval",
  fun = function() {
    build_string_ppi(
      meta_deg$symbol
    )
  },
  attempts = 4L,
  wait_seconds = c(
    10,
    30,
    60
  )
)

mcode_modules <- list(
  `1` = c(
    "PRC1","CENPA","MCM4","MELK","CDKN3","CCNB1","CDK1","CKS2","GINS2","KIF20A",
    "PBK","PTTG1","NCAPG","KIF4A","ATAD2","DLGAP5","CCNA2","UHRF1","RAD51AP1","NUSAP1",
    "CEP55","AURKA","TTK","TPX2","CENPN","TRIP13","ANLN","CDC20","BUB1","CENPW",
    "TOP2A","MAD2L1","ASPM","CDCA7","HJURP","HMMR","FANCI","CDCA5","CENPF","UBE2T","MCM2"
  ),
  `2` = c("CXCL8","SPP1","MMP1","PLAU","COL1A1","CXCL12","MMP3","MMP7","CXCL1","MMP12"),
  `3` = c("MT1H","MT1E","MT1M","MT1G","MT1X","MT2A","MT1F"),
  `4` = c("AQP8","GCG","SLC26A3","GUCA2A","CLCA1","SOX9","CLCA4","MS4A12","GUCA2B"),
  `5` = c("AKR1B10","BCHE","NPY","ADH1B","ADH1A","PYY","VIP","ADH1C","SI","MAOA",
          "GPT","SCG2","GPX3","FABP1","XDH","SST","CHGA","VIPR1")
)

mcode_table <- bind_rows(lapply(names(mcode_modules), function(m) {
  tibble(module = m, symbol = mcode_modules[[m]])
}))
write_tsv_safe(mcode_table, "PPI_MCODE_modules.tsv")

ppi_centrality <- NULL
ppi_metric_cols <- character(0)

# Descending centrality rank + quartile annotation.
# Convention used throughout this analysis:
#   Q1 = highest/most-central 25% of genes
#   Q2 = 25-50%
#   Q3 = 50-75%
#   Q4 = lowest/least-central 25%
# Higher values are treated as stronger centrality for every included metric and
# for the final composite score. Ties receive the same rank/quartile whenever
# possible, so quartile sizes can differ slightly at a boundary.
add_rank_quartile <- function(df, column) {
  x <- as.numeric(df[[column]])
  ok <- is.finite(x)
  
  rank_desc <- rep(NA_integer_, length(x))
  quartile <- rep(NA_character_, length(x))
  q1_or_q2 <- rep(NA, length(x))
  
  if (any(ok)) {
    r <- rank(-x[ok], ties.method = "min")
    n_ok <- sum(ok)
    pct_position <- if (n_ok <= 1) rep(0, n_ok) else (r - 1) / (n_ok - 1)
    
    q <- ifelse(
      pct_position < 0.25, "Q1",
      ifelse(pct_position < 0.50, "Q2",
             ifelse(pct_position < 0.75, "Q3", "Q4"))
    )
    
    rank_desc[ok] <- as.integer(r)
    quartile[ok] <- q
    q1_or_q2[ok] <- q %in% c("Q1", "Q2")
  }
  
  df[[paste0(column, "_rank")]] <- rank_desc
  df[[paste0(column, "_quartile")]] <- quartile
  df[[paste0(column, "_Q1_or_Q2")]] <- q1_or_q2
  df
}

if (!is.null(ppi)) {
  g <- ppi$graph
  write_tsv_safe(ppi$edges, "PPI_STRING_edges.tsv")
  
  exact_cent_file <- file.path(DIR_DATA, "PPI_centrality_export.tsv")
  if (file.exists(exact_cent_file)) {
    ppi_centrality <- readr::read_tsv(exact_cent_file, show_col_types = FALSE)
    if (!"symbol" %in% colnames(ppi_centrality)) {
      stop("Data/PPI_centrality_export.tsv must contain a 'symbol' column.")
    }
  } else {
    # Fully reproducible R fallback. This uses five standard graph metrics.
    # For exact reproduction of the seven-metric CytoNCA/cytoHubba ranking,
    # export that Cytoscape table to Data/PPI_centrality_export.tsv.
    ppi_centrality <- tibble(
      symbol = igraph::V(g)$name,
      degree = igraph::degree(g),
      closeness = igraph::closeness(g, normalized = TRUE),
      betweenness = igraph::betweenness(g, normalized = TRUE),
      eigenvector = igraph::eigen_centrality(g, directed = FALSE, scale = TRUE)$vector,
      clustering_coefficient = igraph::transitivity(g, type = "localundirected", isolates = "zero")
    )
  }
  
  # Identify only genuine numeric centrality metrics. Exclude helper/rank columns
  # if the supplied Cytoscape file already contains them.
  helper_cols <- c("symbol", "name", "module", "score", "final_score", "top20pct")
  candidate_metric_cols <- setdiff(colnames(ppi_centrality), helper_cols)
  candidate_metric_cols <- candidate_metric_cols[
    !grepl("(_rank|_quartile|_Q1_or_Q2)$", candidate_metric_cols)
  ]
  ppi_metric_cols <- candidate_metric_cols[
    vapply(ppi_centrality[candidate_metric_cols], is.numeric, logical(1))
  ]
  
  if (length(ppi_metric_cols) == 0) {
    stop("No numeric PPI centrality metrics were available for ranking.")
  }
  
  # Final composite score: min-max scale every metric to [0,1], protect zeros,
  # multiply across metrics, then log10. Larger final_score = stronger centrality.
  scaled_metrics <- as.data.frame(lapply(ppi_centrality[ppi_metric_cols], minmax01))
  scaled_metrics[] <- lapply(scaled_metrics, function(z) {
    z[!is.finite(z)] <- 0
    z
  })
  composite_product <- apply(pmax(as.matrix(scaled_metrics), 1e-6), 1, prod)
  ppi_centrality$final_score <- log10(composite_product)
  
  # Add rank/quartile information for EACH individual metric.
  for (m in ppi_metric_cols) {
    ppi_centrality <- add_rank_quartile(ppi_centrality, m)
  }
  
  # Add rank/quartile information for the FINAL composite centrality score.
  ppi_centrality <- add_rank_quartile(ppi_centrality, "final_score")
  
  # Compact summary across the individual metrics (final score excluded here).
  metric_quartile_cols <- paste0(ppi_metric_cols, "_quartile")
  ppi_centrality$n_metrics_Q1 <- rowSums(
    ppi_centrality[metric_quartile_cols] == "Q1", na.rm = TRUE
  )
  qmat <- as.matrix(ppi_centrality[metric_quartile_cols])
  ppi_centrality$n_metrics_Q1_or_Q2 <- rowSums(
    qmat == "Q1" | qmat == "Q2", na.rm = TRUE
  )
  
  ppi_centrality <- ppi_centrality |>
    arrange(final_score_rank, symbol)
  
  # Reorder columns so every metric is immediately followed by its rank/quartile.
  metric_blocks <- unlist(lapply(ppi_metric_cols, function(m) {
    c(m, paste0(m, "_rank"), paste0(m, "_quartile"), paste0(m, "_Q1_or_Q2"))
  }))
  ordered_cols <- c(
    "symbol",
    metric_blocks,
    "final_score", "final_score_rank", "final_score_quartile", "final_score_Q1_or_Q2",
    "n_metrics_Q1", "n_metrics_Q1_or_Q2"
  )
  other_cols <- setdiff(colnames(ppi_centrality), ordered_cols)
  ppi_centrality <- ppi_centrality |>
    dplyr::select(any_of(c(ordered_cols, other_cols)))
  
  # Main detailed table and an easy-to-read quartile-only companion table.
  write_tsv_safe(ppi_centrality, "PPI_centrality_rank.tsv")
  
  quartile_only_cols <- c(
    "symbol",
    unlist(lapply(ppi_metric_cols, function(m) {
      c(paste0(m, "_rank"), paste0(m, "_quartile"), paste0(m, "_Q1_or_Q2"))
    })),
    "final_score_rank", "final_score_quartile", "final_score_Q1_or_Q2",
    "n_metrics_Q1", "n_metrics_Q1_or_Q2"
  )
  write_tsv_safe(
    ppi_centrality |> dplyr::select(any_of(quartile_only_cols)),
    "PPI_centrality_quartiles.tsv"
  )
  
  # Publication-quality MCODE module networks are rendered in Section 28.
}

if (!is.null(ppi_centrality) && nrow(ppi_centrality) > 0) {
  
  ppi_top <- ppi_centrality |>
    filter(
      is.finite(final_score),
      !is.na(final_score_rank)
    ) |>
    arrange(final_score_rank) |>
    slice_head(n = 30) |>
    mutate(
      focal = symbol %in% FOCAL_GENES
    )
  
  p_ppi_centrality <- ggplot(
    ppi_top,
    aes(
      x = final_score,
      y = reorder(symbol, final_score),
      color = focal
    )
  ) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c(
        `FALSE` = unname(PLOT_COLORS["down"]),
        `TRUE` = unname(PLOT_COLORS["up"])
      )
    ) +
    labs(
      title = "Highest-ranked STRING/PPI genes",
      x = "Final Centrality score",
      y = NULL,
      color = "Focal gene"
    ) +
    theme(
      legend.position = "bottom"
    )
  
  save_supplementary_plot(
    p_ppi_centrality,
    "Figure_S09_PPI_centrality",
    width = 7.5,
    height = 8
  )
}

###############################################################################
# 21. TCGA DATA DOWNLOAD / CACHE
###############################################################################

download_tcga_project <- function(project) {
  cache_file <- file.path(
    DIR_TCGA,
    paste0(project, "_RNAseq_SE.rds")
  )
  
  if (file.exists(cache_file)) {
    message("Loading cached ", project, " RNA-seq object ...")
    return(readRDS(cache_file))
  }
  
  message("Preparing ", project, " RNA-seq download ...")
  
  q <- retry_remote_call(
    label = paste0(project, " GDC query"),
    fun = function() {
      TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        experimental.strategy = "RNA-Seq",
        workflow.type = "STAR - Counts",
        sample.type = c(
          "Primary Tumor",
          "Solid Tissue Normal"
        )
      )
    },
    attempts = 5L,
    wait_seconds = c(
      15,
      30,
      60,
      120
    )
  )
  
  retry_remote_call(
    label = paste0(project, " GDC download"),
    fun = function() {
      TCGAbiolinks::GDCdownload(
        q,
        method = "api",
        directory = DIR_TCGA,
        files.per.chunk = 40
      )
      TRUE
    },
    attempts = 4L,
    wait_seconds = c(
      20,
      60,
      120
    )
  )
  
  se <- TCGAbiolinks::GDCprepare(
    q,
    directory = DIR_TCGA,
    summarizedExperiment = TRUE
  )
  
  saveRDS(
    se,
    cache_file
  )
  
  se
}

se_coad <- download_tcga_project("TCGA-COAD")
se_read <- download_tcga_project("TCGA-READ")

extract_unstranded <- function(se) {
  an <- SummarizedExperiment::assayNames(se)
  pick <- if ("unstranded" %in% an) "unstranded" else an[1]
  SummarizedExperiment::assay(se, pick)
}

counts_coad <- extract_unstranded(se_coad)
counts_read <- extract_unstranded(se_read)
common_ens <- intersect(rownames(counts_coad), rownames(counts_read))
counts_coad <- counts_coad[common_ens, , drop = FALSE]
counts_read <- counts_read[common_ens, , drop = FALSE]
counts_tcga <- cbind(counts_coad, counts_read)

sample_type_code <- substr(colnames(counts_tcga), 14, 15)
tcga_meta <- tibble(
  sample = colnames(counts_tcga),
  project = c(rep("COAD", ncol(counts_coad)), rep("READ", ncol(counts_read))),
  condition = ifelse(sample_type_code == "11", "normal", "tumor")
) |>
  mutate(
    project = factor(project),
    condition = factor(condition, levels = c("normal", "tumor"))
  )


# -------------------------------------------------------------------------
# TCGA SAMPLE / ALIQUOT AUDIT
# This is reporting-only and does not change the DESeq2 or TPM analyses.
# It records the exact RNA-seq columns downloaded, unique cases represented,
# and any cases with more than one aliquot within the same tissue condition.
# -------------------------------------------------------------------------

tcga_expression_sample_audit <- tcga_meta |>
  dplyr::transmute(
    sample = as.character(sample),
    case_id = substr(as.character(sample), 1, 12),
    project = as.character(project),
    sample_type_code = substr(as.character(sample), 14, 15),
    condition = as.character(condition)
  )

write_tsv_safe(
  tcga_expression_sample_audit,
  "TCGA_expression_sample_audit.tsv"
)

tcga_expression_count_summary <- dplyr::bind_rows(
  tcga_expression_sample_audit |>
    dplyr::count(
      project,
      condition,
      name = "n"
    ) |>
    dplyr::mutate(
      unit = "RNA-seq sample/aliquot columns"
    ),
  tcga_expression_sample_audit |>
    dplyr::distinct(
      case_id,
      project,
      condition
    ) |>
    dplyr::count(
      project,
      condition,
      name = "n"
    ) |>
    dplyr::mutate(
      unit = "unique cases"
    )
) |>
  dplyr::select(
    unit,
    project,
    condition,
    n
  )

write_tsv_safe(
  tcga_expression_count_summary,
  "TCGA_expression_count_summary.tsv"
)

tcga_duplicate_aliquot_audit <- tcga_expression_sample_audit |>
  dplyr::count(
    case_id,
    project,
    condition,
    name = "n_aliquots"
  ) |>
  dplyr::filter(
    n_aliquots > 1
  ) |>
  dplyr::arrange(
    project,
    condition,
    dplyr::desc(n_aliquots),
    case_id
  )

write_tsv_safe(
  tcga_duplicate_aliquot_audit,
  "TCGA_duplicate_aliquot_audit.tsv"
)

tcga_case_tissue_audit <- tcga_expression_sample_audit |>
  dplyr::distinct(
    case_id,
    project,
    condition
  ) |>
  dplyr::mutate(
    present = TRUE
  ) |>
  tidyr::pivot_wider(
    names_from = condition,
    values_from = present,
    values_fill = FALSE
  ) |>
  dplyr::mutate(
    has_tumor = tumor,
    has_normal = normal,
    has_both = tumor & normal
  )

write_tsv_safe(
  tcga_case_tissue_audit,
  "TCGA_case_tissue_audit.tsv"
)

message(
  "TCGA RNA-seq audit: ",
  sum(tcga_expression_sample_audit$condition == "tumor"),
  " tumor columns and ",
  sum(tcga_expression_sample_audit$condition == "normal"),
  " normal columns; ",
  dplyr::n_distinct(
    tcga_expression_sample_audit$case_id[
      tcga_expression_sample_audit$condition == "tumor"
    ]
  ),
  " unique tumor cases and ",
  dplyr::n_distinct(
    tcga_expression_sample_audit$case_id[
      tcga_expression_sample_audit$condition == "normal"
    ]
  ),
  " unique normal cases."
)

# Gene metadata from COAD object; row order restricted to common genes.
gene_info <- as.data.frame(SummarizedExperiment::rowData(se_coad))[common_ens, , drop = FALSE]
gene_symbol <- if ("gene_name" %in% colnames(gene_info)) gene_info$gene_name else rownames(gene_info)
gene_type <- if ("gene_type" %in% colnames(gene_info)) gene_info$gene_type else NA_character_

###############################################################################
# 22. TCGA DIFFERENTIAL EXPRESSION — COAD+READ VALIDATION
###############################################################################

# Collapse multiple aliquots from the same case and tissue condition before
# DESeq2 so they are not treated as independent biological observations.
# Raw counts are summed within case-condition, analogous to collapsing
# technical replicates. The patient-level TPM analysis already averages
# duplicate aliquots separately.

tcga_deseq_group_map <- tcga_expression_sample_audit |>
  dplyr::mutate(
    group_id = paste(
      project,
      case_id,
      condition,
      sep = "__"
    )
  )

counts_tcga_case <- t(
  rowsum(
    t(counts_tcga),
    group = tcga_deseq_group_map$group_id,
    reorder = FALSE
  )
)

tcga_meta_case <- tcga_deseq_group_map |>
  dplyr::distinct(
    group_id,
    .keep_all = TRUE
  ) |>
  dplyr::slice(
    match(
      colnames(counts_tcga_case),
      group_id
    )
  ) |>
  dplyr::mutate(
    project = factor(project),
    condition = factor(
      condition,
      levels = c(
        "normal",
        "tumor"
      )
    )
  )

if (
  anyNA(tcga_meta_case$group_id) ||
  !identical(
    colnames(counts_tcga_case),
    tcga_meta_case$group_id
  )
) {
  stop(
    "TCGA case-level count collapsing produced a metadata alignment error."
  )
}

tcga_deseq_case_summary <- tcga_meta_case |>
  dplyr::count(
    project,
    condition,
    name = "n_case_condition_profiles"
  )

write_tsv_safe(
  tcga_deseq_case_summary,
  "TCGA_DESeq2_case_level_input_summary.tsv"
)

write_tsv_safe(
  tcga_meta_case |>
    dplyr::select(
      group_id,
      case_id,
      project,
      condition,
      sample_type_code
    ),
  "TCGA_DESeq2_case_level_input_metadata.tsv"
)

# -------------------------------------------------------------------------
# Memory-safe DESeq2 preparation
#
# Important:
#   - The statistical model is unchanged: ~ project + condition.
#   - The original low-count rule is unchanged: total count >= 10.
#   - Filtering is performed BEFORE constructing the DESeqDataSet to avoid
#     carrying tens of thousands of uninformative rows through DESeq2.
#   - Large SummarizedExperiment/count objects are temporarily released
#     before DESeq2 and reloaded from the local cache afterwards.
# -------------------------------------------------------------------------

tcga_gene_keep <- rowSums(
  counts_tcga_case,
  na.rm = TRUE
) >= 10

message(
  "TCGA DESeq2 prefilter: retaining ",
  sum(tcga_gene_keep),
  " / ",
  length(tcga_gene_keep),
  " genes with total count >= 10."
)

write_tsv_safe(
  tibble::tibble(
    ensembl = rownames(counts_tcga_case),
    total_count = rowSums(
      counts_tcga_case,
      na.rm = TRUE
    ),
    retained_for_DESeq2 = tcga_gene_keep
  ),
  "TCGA_DESeq2_prefilter_audit.tsv"
)

counts_tcga_case <- counts_tcga_case[
  tcga_gene_keep,
  ,
  drop = FALSE
]

# Counts are sums of integer STAR counts. Round once, then store as integer so
# DESeqDataSetFromMatrix does not need to create another full rounded copy.
counts_tcga_case <- round(
  counts_tcga_case
)

if (
  anyNA(counts_tcga_case) ||
  any(counts_tcga_case < 0)
) {
  stop(
    "TCGA case-level count matrix contains missing or negative values."
  )
}

storage.mode(counts_tcga_case) <- "integer"

# Release the much larger source objects while DESeq2 is fitting.
# The cached SummarizedExperiment objects are reloaded below for TPM/clinical
# analyses, so this does not change the data used later in the pipeline.
rm(
  counts_coad,
  counts_read,
  counts_tcga,
  se_coad,
  se_read
)
invisible(
  gc(
    full = TRUE
  )
)

dds_tcga <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_tcga_case,
  colData = data.frame(
    row.names = tcga_meta_case$group_id,
    project = tcga_meta_case$project,
    condition = tcga_meta_case$condition
  ),
  design = ~ project + condition
)

# Once the DESeqDataSet owns the count matrix, drop the standalone copy.
rm(
  counts_tcga_case
)
invisible(
  gc(
    full = TRUE
  )
)

dds_tcga <- tryCatch(
  DESeq2::DESeq(
    dds_tcga,
    quiet = TRUE,
    parallel = FALSE
  ),
  error = function(e) {
    msg <- conditionMessage(e)

    if (
      grepl(
        "bad_alloc|cannot allocate|memory",
        msg,
        ignore.case = TRUE
      )
    ) {
      stop(
        paste0(
          "DESeq2 stopped because R could not allocate enough memory even ",
          "after pre-filtering and releasing the cached TCGA objects. ",
          "Close other R sessions/applications and rerun from the TCGA section. ",
          "Original error: ",
          msg
        ),
        call. = FALSE
      )
    }

    stop(e)
  }
)

res_tcga <- DESeq2::results(
  dds_tcga,
  contrast = c(
    "condition",
    "tumor",
    "normal"
  )
)

tcga_de <- as.data.frame(res_tcga) |>
  rownames_to_column("ensembl") |>
  mutate(
    symbol = gene_symbol[match(ensembl, common_ens)],
    gene_type = gene_type[match(ensembl, common_ens)]
  ) |>
  arrange(padj)

write_tsv_safe(tcga_de, "TCGA_COAD_READ_DESeq2_tumor_vs_normal.tsv")

# Combined TCGA-COAD + TCGA-READ top-25 genes by tumor-vs-normal effect size.
# This uses the same joint DESeq2 model above (~ project + condition), so the ranking
# reflects the pooled COAD+READ tumor effect while adjusting for project.
tcga_coad_read_ranked <- tcga_de |>
  filter(!is.na(symbol), !is.na(padj), is.finite(log2FoldChange))

tcga_coad_read_top25_up <- tcga_coad_read_ranked |>
  filter(padj < 0.05) |>
  arrange(desc(log2FoldChange), padj) |>
  distinct(symbol, .keep_all = TRUE) |>
  slice_head(n = 25)

tcga_coad_read_top25_down <- tcga_coad_read_ranked |>
  filter(padj < 0.05) |>
  arrange(log2FoldChange, padj) |>
  distinct(symbol, .keep_all = TRUE) |>
  slice_head(n = 25)

write_tsv_safe(tcga_coad_read_top25_up, "TCGA_COAD_READ_top25_up.tsv")
write_tsv_safe(tcga_coad_read_top25_down, "TCGA_COAD_READ_top25_down.tsv")

tcga_volcano_df <- tcga_de |>
  filter(
    is.finite(log2FoldChange),
    !is.na(padj),
    is.finite(padj)
  ) |>
  mutate(
    minus_log10_fdr = -log10(
      pmax(
        padj,
        .Machine$double.xmin
      )
    ),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "Not significant"
    ),
    label = ifelse(
      symbol %in% FOCAL_GENES,
      symbol,
      NA_character_
    )
  )

p_tcga_volcano <- ggplot(
  tcga_volcano_df,
  aes(
    x = log2FoldChange,
    y = minus_log10_fdr,
    color = direction
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggrepel::geom_label_repel(
    data = subset(
      tcga_volcano_df,
      !is.na(label)
    ),
    aes(label = label),
    show.legend = FALSE,
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      Down = unname(PLOT_COLORS["down"]),
      `Not significant` = unname(PLOT_COLORS["neutral"]),
      Up = unname(PLOT_COLORS["up"])
    )
  ) +
  labs(
    title = "TCGA-COAD/READ differential expression",
    x = expression("log"[2]*" fold change (tumor vs normal)"),
    y = expression(-log[10](FDR)),
    color = NULL
  ) +
  theme(
    legend.position = "bottom"
  )

tcga_top25_plot_df <- bind_rows(
  tcga_coad_read_top25_up |>
    mutate(direction = "Up"),
  tcga_coad_read_top25_down |>
    mutate(direction = "Down")
) |>
  mutate(
    direction = factor(
      direction,
      levels = c("Down", "Up")
    )
  )

p_tcga_top25 <- ggplot(
  tcga_top25_plot_df,
  aes(
    x = log2FoldChange,
    y = reorder(symbol, log2FoldChange),
    fill = direction
  )
) +
  geom_col(width = 0.75) +
  facet_wrap(
    ~direction,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      Down = unname(PLOT_COLORS["down"]),
      Up = unname(PLOT_COLORS["up"])
    )
  ) +
  labs(
    title = "Top TCGA tumor-versus-normal expression changes",
    x = expression("log"[2]*" fold change"),
    y = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 7.5)
  )

save_supplementary_plot(
  p_tcga_volcano /
    p_tcga_top25 +
    patchwork::plot_annotation(tag_levels = "A"),
  "Figure_S10_TCGA_differential_expression",
  width = 12,
  height = 12
)


# DESeq2 outputs needed downstream have now been materialized as ordinary
# tables/plots. Release the fitted DESeqDataSet before reloading TPM assays.
if (exists("dds_tcga")) {
  rm(dds_tcga)
}
if (exists("res_tcga")) {
  rm(res_tcga)
}
invisible(
  gc(
    full = TRUE
  )
)


###############################################################################
# 23. TCGA FOCAL-GENE EXPRESSION
###############################################################################

# Focal-gene validation is based on log2(TPM + 1), using the same transcript
# selection and patient-level aggregation as the survival analysis. The TPM
# matrix and publication plots are constructed below after clinical preparation.


###############################################################################
# 24. TCGA SURVIVAL ANALYSIS
###############################################################################

# Reload the cached TCGA objects that were deliberately released before DESeq2
# to reduce peak memory use. This is a local cache read, not a new download.
if (!exists("se_coad")) {
  se_coad <- download_tcga_project("TCGA-COAD")
}

if (!exists("se_read")) {
  se_read <- download_tcga_project("TCGA-READ")
}

invisible(
  gc(
    full = TRUE
  )
)

as_numeric_safe <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

is_missing_text <- function(x) {
  y <- trimws(
    as.character(x)
  )
  
  is.na(y) |
    y == "" |
    tolower(y) %in% c(
      "na",
      "n/a",
      "not reported",
      "not available",
      "unknown",
      "[not available]",
      "[not reported]",
      "--"
    )
}

candidate_columns <- function(
    df,
    exact = character(0),
    regex = character(0)
) {
  
  out <- exact[
    exact %in% colnames(df)
  ]
  
  if (length(regex) > 0) {
    for (pattern in regex) {
      out <- c(
        out,
        grep(
          pattern,
          colnames(df),
          ignore.case = TRUE,
          value = TRUE
        )
      )
    }
  }
  
  unique(out)
}

coalesce_character <- function(
    df,
    exact = character(0),
    regex = character(0)
) {
  
  columns <- candidate_columns(
    df,
    exact,
    regex
  )
  
  out <- rep(
    NA_character_,
    nrow(df)
  )
  
  for (column in columns) {
    
    x <- as.character(
      df[[column]]
    )
    
    x[is_missing_text(x)] <- NA_character_
    
    use <- is.na(out) & !is.na(x)
    
    out[use] <- x[use]
  }
  
  out
}

coalesce_numeric <- function(
    df,
    exact = character(0),
    regex = character(0)
) {
  
  columns <- candidate_columns(
    df,
    exact,
    regex
  )
  
  out <- rep(
    NA_real_,
    nrow(df)
  )
  
  for (column in columns) {
    
    x <- as_numeric_safe(
      df[[column]]
    )
    
    use <- !is.finite(out) & is.finite(x)
    
    out[use] <- x[use]
  }
  
  out
}

first_nonmissing_character <- function(x) {
  
  x <- as.character(x)
  x[is_missing_text(x)] <- NA_character_
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    NA_character_
  } else {
    x[1]
  }
}

first_finite_numeric <- function(x) {
  
  x <- as_numeric_safe(x)
  x <- x[is.finite(x)]
  
  if (length(x) == 0) {
    NA_real_
  } else {
    x[1]
  }
}

normalize_pathologic_stage <- function(x) {
  
  y <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  y[is_missing_text(y)] <- NA_character_
  
  y <- gsub(
    "^AJCC[ _-]*PATHOLOGIC(AL)?[ _-]*STAGE[ _-]*",
    "",
    y
  )
  
  y <- gsub(
    "^PATHOLOGIC(AL)?[ _-]*STAGE[ _-]*",
    "",
    y
  )
  
  y <- gsub(
    "^STAGE[ _-]*",
    "",
    y
  )
  
  y <- gsub(
    "[^A-Z0-9]",
    "",
    y
  )
  
  stage <- case_when(
    grepl("^IV", y) ~ "IV",
    grepl("^III", y) ~ "III",
    grepl("^II", y) ~ "II",
    grepl("^I", y) ~ "I",
    TRUE ~ NA_character_
  )
  
  factor(
    stage,
    levels = c(
      "I",
      "II",
      "III",
      "IV"
    )
  )
}

clinical_from_se <- function(
    se,
    project_label
) {
  
  cd <- as.data.frame(
    SummarizedExperiment::colData(se),
    stringsAsFactors = FALSE
  )
  
  cd$sample <- rownames(cd)
  
  vital_status <- coalesce_character(
    cd,
    exact = c(
      "vital_status",
      "paper_vital_status"
    ),
    regex = c(
      "(^|[._])vital_status$"
    )
  )
  
  days_to_death <- coalesce_numeric(
    cd,
    exact = c(
      "days_to_death",
      "paper_days_to_death"
    ),
    regex = c(
      "(^|[._])days_to_death$"
    )
  )
  
  days_to_last_follow_up <- coalesce_numeric(
    cd,
    exact = c(
      "days_to_last_follow_up",
      "days_to_last_followup",
      "days_to_last_known_alive",
      "paper_days_to_last_follow_up"
    ),
    regex = c(
      "days_to_last_follow",
      "days_to_last_known_alive"
    )
  )
  
  age_at_index <- coalesce_numeric(
    cd,
    exact = c(
      "age_at_index",
      "paper_age_at_index"
    ),
    regex = c(
      "(^|[._])age_at_index$"
    )
  )
  
  age_at_diagnosis <- coalesce_numeric(
    cd,
    exact = c(
      "age_at_diagnosis",
      "paper_age_at_diagnosis"
    ),
    regex = c(
      "(^|[._])age_at_diagnosis$"
    )
  )
  
  gender <- coalesce_character(
    cd,
    exact = c(
      "gender",
      "sex",
      "paper_gender"
    ),
    regex = c(
      "(^|[._])gender$",
      "(^|[._])sex$"
    )
  )
  
  pathologic_stage <- coalesce_character(
    cd,
    exact = c(
      "ajcc_pathologic_stage",
      "pathologic_stage",
      "pathological_stage",
      "tumor_stage",
      "paper_pathologic_stage",
      "paper_stage_event_pathologic_stage"
    ),
    regex = c(
      "ajcc.*pathologic.*stage",
      "pathologic.*stage",
      "pathological.*stage",
      "(^|[._])tumor_stage$"
    )
  )
  
  tibble(
    sample = cd$sample,
    case_id = substr(
      cd$sample,
      1,
      12
    ),
    project = project_label,
    vital_status = vital_status,
    days_to_death = days_to_death,
    days_to_last_follow_up = days_to_last_follow_up,
    age_at_index = age_at_index,
    age_at_diagnosis = age_at_diagnosis,
    sex = gender,
    raw_stage = pathologic_stage
  ) |>
    group_by(
      case_id,
      project
    ) |>
    summarise(
      vital_status = first_nonmissing_character(vital_status),
      days_to_death = first_finite_numeric(days_to_death),
      days_to_last_follow_up = first_finite_numeric(days_to_last_follow_up),
      age_at_index = first_finite_numeric(age_at_index),
      age_at_diagnosis = first_finite_numeric(age_at_diagnosis),
      sex = first_nonmissing_character(sex),
      raw_stage = first_nonmissing_character(raw_stage),
      .groups = "drop"
    )
}

clinical_tcga <- bind_rows(
  clinical_from_se(
    se_coad,
    "COAD"
  ),
  clinical_from_se(
    se_read,
    "READ"
  )
)

write_tsv_safe(
  clinical_tcga,
  "TCGA_clinical_from_SummarizedExperiment.tsv"
)


# -------------------------------------------------------------------------
# TCGA SEX / GENDER METADATA AUDIT
#
# Two sources are compared:
#   1) the metadata embedded in the downloaded SummarizedExperiment objects;
#   2) a fresh case-level GDC clinical query through TCGAbiolinks.
#
# The final survival models use case-level GDC sex_at_birth as an analytical
# covariate. This is therefore a required live input, not a reporting-only
# audit. If the GDC clinical endpoint cannot be queried, the pipeline stops
# rather than falling back to incomplete embedded sex/gender metadata.
# -------------------------------------------------------------------------

standardize_sex_audit <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y[is_missing_text(y)] <- NA_character_
  dplyr::case_when(
    y %in% c("male", "m", "man") ~ "male",
    y %in% c("female", "f", "woman") ~ "female",
    TRUE ~ y
  )
}

sex_se_audit <- clinical_tcga |>
  dplyr::transmute(
    case_id = as.character(case_id),
    project = as.character(project),
    sex_se_raw = as.character(sex),
    sex_se = standardize_sex_audit(sex)
  )

write_tsv_safe(
  sex_se_audit,
  "TCGA_sex_SummarizedExperiment_audit.tsv"
)

sex_se_summary <- sex_se_audit |>
  dplyr::group_by(project) |>
  dplyr::summarise(
    n_cases = dplyr::n_distinct(case_id),
    n_sex_available = sum(!is.na(sex_se)),
    proportion_available = n_sex_available / n_cases,
    .groups = "drop"
  )

write_tsv_safe(
  sex_se_summary,
  "TCGA_sex_SummarizedExperiment_summary.tsv"
)

extract_gdc_sex_value <- function(demographic) {

  if (is.null(demographic)) {
    return(NA_character_)
  }

  if (
    is.list(demographic) &&
    !is.null(demographic$sex_at_birth)
  ) {
    return(
      as.character(
        demographic$sex_at_birth
      )[1]
    )
  }

  if (
    is.list(demographic) &&
    length(demographic) > 0 &&
    is.list(demographic[[1]]) &&
    !is.null(demographic[[1]]$sex_at_birth)
  ) {
    return(
      as.character(
        demographic[[1]]$sex_at_birth
      )[1]
    )
  }

  NA_character_
}

query_gdc_sex_audit <- function(
    project_full,
    project_short
) {

  gdc_filter <- list(
    op = "in",
    content = list(
      field = "cases.project.project_id",
      value = list(
        project_full
      )
    )
  )

  response <- httr2::request(
    "https://api.gdc.cancer.gov/cases"
  ) |>
    httr2::req_url_query(
      filters = jsonlite::toJSON(
        gdc_filter,
        auto_unbox = TRUE
      ),
      fields = paste(
        c(
          "submitter_id",
          "project.project_id",
          "demographic.sex_at_birth"
        ),
        collapse = ","
      ),
      expand = "demographic",
      format = "JSON",
      size = 2000
    ) |>
    httr2::req_retry(
      max_tries = 4
    ) |>
    httr2::req_timeout(
      seconds = 120
    ) |>
    httr2::req_perform()

  body <- httr2::resp_body_json(
    response,
    simplifyVector = FALSE
  )

  hits <- body$data$hits

  if (
    is.null(hits) ||
    length(hits) == 0
  ) {
    stop(
      project_full,
      ": the GDC Cases API returned no cases."
    )
  }

  dplyr::bind_rows(
    lapply(
      hits,
      function(hit) {
        tibble::tibble(
          case_id = substr(
            as.character(
              hit$submitter_id
            ),
            1,
            12
          ),
          project = project_short,
          sex_gdc_raw = extract_gdc_sex_value(
            hit$demographic
          )
        )
      }
    )
  ) |>
    dplyr::mutate(
      sex_gdc = standardize_sex_audit(
        sex_gdc_raw
      ),
      source_case_column = "submitter_id",
      source_sex_column = "demographic.sex_at_birth"
    ) |>
    dplyr::filter(
      grepl(
        "^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}$",
        case_id
      )
    ) |>
    dplyr::group_by(
      case_id,
      project
    ) |>
    dplyr::summarise(
      sex_gdc_raw = first_nonmissing_character(
        sex_gdc_raw
      ),
      sex_gdc = first_nonmissing_character(
        sex_gdc
      ),
      source_case_column = first_nonmissing_character(
        source_case_column
      ),
      source_sex_column = first_nonmissing_character(
        source_sex_column
      ),
      .groups = "drop"
    )
}


gdc_sex_audit <- retry_remote_call(
  label = "GDC case-level sex_at_birth query",
  fun = function() {
    dplyr::bind_rows(
      query_gdc_sex_audit(
        "TCGA-COAD",
        "COAD"
      ),
      query_gdc_sex_audit(
        "TCGA-READ",
        "READ"
      )
    )
  },
  attempts = 4L,
  wait_seconds = c(
    10,
    30,
    60
  )
)

if (!is.null(gdc_sex_audit)) {

  write_tsv_safe(
    gdc_sex_audit,
    "TCGA_sex_GDC_clinical_audit.tsv"
  )

  gdc_sex_summary <- gdc_sex_audit |>
    dplyr::group_by(project) |>
    dplyr::summarise(
      n_cases = dplyr::n_distinct(case_id),
      n_sex_available = sum(!is.na(sex_gdc)),
      proportion_available = n_sex_available / n_cases,
      .groups = "drop"
    )

  write_tsv_safe(
    gdc_sex_summary,
    "TCGA_sex_GDC_clinical_summary.tsv"
  )

  sex_source_comparison <- dplyr::full_join(
    sex_se_audit,
    gdc_sex_audit,
    by = c(
      "case_id",
      "project"
    )
  ) |>
    dplyr::mutate(
      comparison = dplyr::case_when(
        is.na(sex_se) & is.na(sex_gdc) ~
          "missing in both",
        is.na(sex_se) & !is.na(sex_gdc) ~
          "missing in SummarizedExperiment only",
        !is.na(sex_se) & is.na(sex_gdc) ~
          "missing in live GDC only",
        sex_se == sex_gdc ~
          "match",
        TRUE ~
          "discordant"
      )
    ) |>
    dplyr::arrange(
      project,
      case_id
    )

  write_tsv_safe(
    sex_source_comparison,
    "TCGA_sex_source_comparison.tsv"
  )

  sex_source_summary <- sex_source_comparison |>
    dplyr::count(
      project,
      comparison,
      name = "n_cases"
    )

  write_tsv_safe(
    sex_source_summary,
    "TCGA_sex_source_comparison_summary.tsv"
  )

  message(
    "TCGA sex audit — SummarizedExperiment available: ",
    sum(!is.na(sex_source_comparison$sex_se)),
    "/",
    nrow(sex_source_comparison),
    "; live GDC available: ",
    sum(!is.na(sex_source_comparison$sex_gdc)),
    "/",
    nrow(sex_source_comparison),
    "; discordant: ",
    sum(
      sex_source_comparison$comparison == "discordant",
      na.rm = TRUE
    ),
    "."
  )
}


clinical_stage_data <- clinical_tcga |>
  mutate(
    stage_major = normalize_pathologic_stage(
      raw_stage
    ),
    project = factor(
      project
    )
  ) |>
  dplyr::select(
    case_id,
    project,
    stage_major
  ) |>
  filter(
    !is.na(stage_major),
    as.character(stage_major) %in%
      c(
        "I",
        "II",
        "III",
        "IV"
      )
  )

clinical_survival_source <- clinical_tcga

if (
  exists("gdc_sex_audit") &&
  !is.null(gdc_sex_audit)
) {
  clinical_survival_source <- clinical_survival_source |>
    dplyr::left_join(
      gdc_sex_audit |>
        dplyr::select(
          case_id,
          project,
          sex_gdc
        ),
      by = c(
        "case_id",
        "project"
      )
    ) |>
    dplyr::mutate(
      sex = sex_gdc
    ) |>
    dplyr::select(
      -sex_gdc
    )
}

clinical_survival <- clinical_survival_source |>
  mutate(
    vital_status_clean = tolower(
      trimws(
        vital_status
      )
    ),
    
    event = case_when(
      vital_status_clean == "dead" ~ 1L,
      vital_status_clean == "alive" ~ 0L,
      is.finite(days_to_death) &
        days_to_death > 0 ~ 1L,
      TRUE ~ NA_integer_
    ),
    
    OS_days = case_when(
      event == 1L &
        is.finite(days_to_death) &
        days_to_death > 0 ~ days_to_death,
      
      is.finite(days_to_last_follow_up) &
        days_to_last_follow_up > 0 ~ days_to_last_follow_up,
      
      is.finite(days_to_death) &
        days_to_death > 0 ~ days_to_death,
      
      TRUE ~ NA_real_
    ),
    
    age = case_when(
      is.finite(age_at_index) &
        age_at_index > 0 ~ age_at_index,
      
      is.finite(age_at_diagnosis) &
        abs(age_at_diagnosis) > 365 ~
        abs(age_at_diagnosis) / 365.25,
      
      is.finite(age_at_diagnosis) &
        age_at_diagnosis > 0 ~ age_at_diagnosis,
      
      TRUE ~ NA_real_
    ),
    
    sex = factor(
      standardize_sex_audit(sex),
      levels = c(
        "female",
        "male"
      )
    ),
    
    stage_major = normalize_pathologic_stage(
      raw_stage
    ),
    
    project = factor(
      project
    )
  ) |>
  dplyr::select(
    case_id,
    project,
    event,
    OS_days,
    age,
    sex,
    stage_major
  ) |>
  filter(
    !is.na(event),
    is.finite(OS_days),
    OS_days > 0
  )

write_tsv_safe(
  clinical_survival,
  "TCGA_survival_clinical_patient_level.tsv"
)


tcga_survival_sex_model_audit <- clinical_survival |>
  dplyr::count(
    project,
    sex,
    name = "n_patients"
  ) |>
  dplyr::arrange(
    project,
    sex
  )

write_tsv_safe(
  tcga_survival_sex_model_audit,
  "TCGA_survival_sex_model_audit.tsv"
)


# Compare sex availability in the exact survival-eligible cohort against
# the fresh case-level GDC sex-at-birth query.
if (
  exists("gdc_sex_audit") &&
  !is.null(gdc_sex_audit)
) {

  tcga_survival553_sex_comparison <- clinical_survival |>
    dplyr::transmute(
      case_id = as.character(case_id),
      project = as.character(project),
      sex_locked_pipeline = standardize_sex_audit(
        sex
      )
    ) |>
    dplyr::left_join(
      gdc_sex_audit |>
        dplyr::select(
          case_id,
          project,
          sex_gdc
        ),
      by = c(
        "case_id",
        "project"
      )
    ) |>
    dplyr::mutate(
      comparison = dplyr::case_when(
        is.na(sex_locked_pipeline) &
          is.na(sex_gdc) ~
          "missing in both",
        is.na(sex_locked_pipeline) &
          !is.na(sex_gdc) ~
          "missing in locked pipeline only",
        !is.na(sex_locked_pipeline) &
          is.na(sex_gdc) ~
          "missing in live GDC only",
        sex_locked_pipeline == sex_gdc ~
          "match",
        TRUE ~
          "discordant"
      )
    )

  write_tsv_safe(
    tcga_survival553_sex_comparison,
    "TCGA_survival553_sex_comparison.tsv"
  )

  tcga_survival553_sex_summary <- tcga_survival553_sex_comparison |>
    dplyr::summarise(
      n_survival_eligible = dplyr::n(),
      n_sex_locked_pipeline = sum(
        !is.na(sex_locked_pipeline)
      ),
      n_sex_live_GDC = sum(
        !is.na(sex_gdc)
      ),
      n_matches = sum(
        comparison == "match",
        na.rm = TRUE
      ),
      n_discordant = sum(
        comparison == "discordant",
        na.rm = TRUE
      ),
      n_missing_locked_only = sum(
        comparison == "missing in locked pipeline only",
        na.rm = TRUE
      ),
      n_missing_live_GDC_only = sum(
        comparison == "missing in live GDC only",
        na.rm = TRUE
      )
    )

  write_tsv_safe(
    tcga_survival553_sex_summary,
    "TCGA_survival553_sex_summary.tsv"
  )
}

clinical_completeness <- tibble(
  variable = c(
    "event",
    "OS_days",
    "age",
    "sex",
    "stage_major"
  ),
  available = c(
    sum(!is.na(clinical_survival$event)),
    sum(is.finite(clinical_survival$OS_days)),
    sum(is.finite(clinical_survival$age)),
    sum(!is.na(clinical_survival$sex)),
    sum(!is.na(clinical_survival$stage_major))
  ),
  total = nrow(clinical_survival)
) |>
  mutate(
    proportion = available / total
  )

write_tsv_safe(
  clinical_completeness,
  "TCGA_survival_clinical_completeness.tsv"
)


# TPM expression for survival

required_tpm_assay <- "tpm_unstrand"

if (
  !required_tpm_assay %in%
  SummarizedExperiment::assayNames(se_coad) ||
  !required_tpm_assay %in%
  SummarizedExperiment::assayNames(se_read)
) {
  stop(
    "The tpm_unstrand assay is required for the survival analysis."
  )
}

tpm_coad <- SummarizedExperiment::assay(
  se_coad,
  required_tpm_assay
)

tpm_read <- SummarizedExperiment::assay(
  se_read,
  required_tpm_assay
)

common_tpm_genes <- intersect(
  rownames(tpm_coad),
  rownames(tpm_read)
)

tpm_all <- cbind(
  tpm_coad[
    common_tpm_genes,
    ,
    drop = FALSE
  ],
  tpm_read[
    common_tpm_genes,
    ,
    drop = FALSE
  ]
)

tpm_gene_info <- as.data.frame(
  SummarizedExperiment::rowData(se_coad)
)[
  common_tpm_genes,
  ,
  drop = FALSE
]

if (
  !"gene_name" %in%
  colnames(tpm_gene_info)
) {
  stop(
    "gene_name is missing from rowData(se_coad)."
  )
}

tpm_symbols <- as.character(
  tpm_gene_info$gene_name
)

names(tpm_symbols) <- common_tpm_genes

primary_tumor_samples <- colnames(tpm_all)[
  substr(
    colnames(tpm_all),
    14,
    15
  ) == "01"
]

tpm_tumor <- tpm_all[
  ,
  primary_tumor_samples,
  drop = FALSE
]

best_tpm_row <- function(gene) {
  
  candidates <- names(tpm_symbols)[
    tpm_symbols == gene
  ]
  
  candidates <- intersect(
    candidates,
    rownames(tpm_tumor)
  )
  
  if (length(candidates) == 0) {
    stop(
      "No TPM row found for ",
      gene
    )
  }
  
  if (length(candidates) == 1) {
    return(candidates)
  }
  
  mean_tpm <- rowMeans(
    tpm_tumor[
      candidates,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )
  
  candidates[
    which.max(mean_tpm)
  ]
}

tpm_rows <- setNames(
  vapply(
    FOCAL_GENES,
    best_tpm_row,
    character(1)
  ),
  FOCAL_GENES
)

tpm_expression_long <- bind_rows(
  lapply(
    names(tpm_rows),
    function(gene_name) {
      
      row_id <- unname(
        tpm_rows[gene_name]
      )
      
      if (
        length(row_id) != 1L ||
        is.na(row_id) ||
        !row_id %in% rownames(tpm_tumor)
      ) {
        stop(
          "Could not identify a unique TPM row for ",
          gene_name
        )
      }
      
      expression_values <- log2(
        as.numeric(
          tpm_tumor[
            row_id,
            primary_tumor_samples
          ]
        ) + 1
      )
      
      tibble(
        sample = primary_tumor_samples,
        case_id = substr(
          primary_tumor_samples,
          1,
          12
        ),
        gene = gene_name,
        expression = expression_values
      )
    }
  )
)

tpm_expression_patient <- tpm_expression_long |>
  group_by(
    case_id,
    gene
  ) |>
  summarise(
    expression = mean(
      expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = gene,
    values_from = expression
  )

survival_data_tpm <- inner_join(
  clinical_survival,
  tpm_expression_patient,
  by = "case_id"
) |>
  filter(
    !is.na(stage_major),
    as.character(stage_major) %in%
      c(
        "I",
        "II",
        "III",
        "IV"
      )
  ) |>
  mutate(
    stage_group = factor(
      if_else(
        as.character(stage_major) %in%
          c(
            "I",
            "II"
          ),
        "Stage_I_II",
        "Stage_III_IV"
      ),
      levels = c(
        "Stage_I_II",
        "Stage_III_IV"
      )
    )
  )

write_tsv_safe(
  survival_data_tpm,
  "TCGA_COAD_READ_survival_log2TPM_dataset.tsv"
)

write_tsv_safe(
  survival_data_tpm,
  "TCGA_survival_analysis_dataset.tsv"
)

survival_cohort_audit <- survival_data_tpm |>
  group_by(
    stage_group
  ) |>
  summarise(
    n = n(),
    events = sum(
      event,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_tsv_safe(
  survival_cohort_audit,
  "TCGA_TPM_stage_cohort_audit.tsv"
)

message(
  "Stage-known TPM survival cohort: ",
  nrow(survival_data_tpm),
  " patients; events=",
  sum(
    survival_data_tpm$event,
    na.rm = TRUE
  )
)


# Tumor expression by stage group

tcga_stage_expression_data <- dplyr::inner_join(
  clinical_stage_data,
  tpm_expression_patient,
  by = "case_id"
) |>
  dplyr::mutate(
    stage_group = factor(
      dplyr::if_else(
        as.character(
          stage_major
        ) %in%
          c(
            "I",
            "II"
          ),
        "Stage_I_II",
        "Stage_III_IV"
      ),
      levels = c(
        "Stage_I_II",
        "Stage_III_IV"
      )
    )
  )

write_tsv_safe(
  tcga_stage_expression_data,
  "TCGA_focal_gene_expression_by_stage_group_data.tsv"
)


# Overall continuous Cox model

fit_overall_cox <- function(
    gene,
    data
) {
  
  d <- data |>
    filter(
      !is.na(.data[[gene]]),
      is.finite(.data[[gene]]),
      is.finite(age),
      !is.na(sex),
      !is.na(project),
      !is.na(stage_major)
    )
  
  expression_mean <- mean(
    d[[gene]],
    na.rm = TRUE
  )
  
  expression_sd <- sd(
    d[[gene]],
    na.rm = TRUE
  )
  
  if (
    !is.finite(expression_sd) ||
    expression_sd <= 0
  ) {
    return(NULL)
  }
  
  d <- d |>
    mutate(
      expr_z =
        (.data[[gene]] - expression_mean) /
        expression_sd
    )
  
  fit <- survival::coxph(
    survival::Surv(
      OS_days,
      event
    ) ~
      expr_z +
      age +
      sex +
      project +
      stage_major,
    data = d,
    x = TRUE,
    model = TRUE
  )
  
  sm <- summary(fit)
  
  ci <- exp(
    confint(fit)[
      "expr_z",
    ]
  )
  
  ph <- tryCatch(
    survival::cox.zph(
      fit
    ),
    error = function(e) NULL
  )
  
  ph_gene <- NA_real_
  ph_global <- NA_real_
  
  if (!is.null(ph)) {
    
    if (
      "expr_z" %in%
      rownames(ph$table)
    ) {
      ph_gene <- ph$table[
        "expr_z",
        "p"
      ]
    }
    
    if (
      "GLOBAL" %in%
      rownames(ph$table)
    ) {
      ph_global <- ph$table[
        "GLOBAL",
        "p"
      ]
    }
  }
  
  list(
    fit = fit,
    summary = tibble(
      gene = gene,
      n = nrow(d),
      events = sum(
        d$event,
        na.rm = TRUE
      ),
      HR_per_SD =
        sm$coefficients[
          "expr_z",
          "exp(coef)"
        ],
      CI_low = ci[1],
      CI_high = ci[2],
      p_value =
        sm$coefficients[
          "expr_z",
          "Pr(>|z|)"
        ],
      PH_expression_p = ph_gene,
      PH_global_p = ph_global
    )
  )
}

overall_cox_fits <- lapply(
  FOCAL_GENES,
  fit_overall_cox,
  data = survival_data_tpm
)

names(overall_cox_fits) <- FOCAL_GENES

overall_cox_fits <- overall_cox_fits[
  !vapply(
    overall_cox_fits,
    is.null,
    logical(1)
  )
]

overall_cox_summary <- bind_rows(
  lapply(
    overall_cox_fits,
    `[[`,
    "summary"
  )
) |>
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )

write_tsv_safe(
  overall_cox_summary,
  "TCGA_TPM_overall_continuous_Cox_focal_genes.tsv"
)

p_overall_cox <- ggplot(
  overall_cox_summary,
  aes(
    x = HR_per_SD,
    y = reorder(
      gene,
      HR_per_SD
    ),
    color = gene
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  scale_color_manual(
    values = GENE_COLORS
  ) +
  scale_x_log10() +
  labs(
    title = "Overall adjusted survival associations",
    subtitle = "Continuous tumor expression; age, sex, project, and pathological stage adjusted",
    x = "Hazard ratio per 1-SD increase in log2(TPM+1)",
    y = NULL,
    color = NULL
  ) +
  theme(
    legend.position = "none"
  )

save_plot(
  p_overall_cox,
  "TCGA_TPM_overall_Cox_forest",
  width = 7,
  height = 4.5
)

save_supplementary_plot(
  p_overall_cox,
  "Figure_S11_TCGA_overall_Cox",
  width = 7,
  height = 4.5
)


# Stage-dependent Cox model
# Reduced: Surv ~ expr_z + age + sex + project + stage_major
# Full:    Surv ~ expr_z + expr_z:stage_group + age + sex + project + stage_major

fit_stage_interaction <- function(
    gene,
    data
) {
  
  d <- data |>
    filter(
      !is.na(.data[[gene]]),
      is.finite(.data[[gene]]),
      is.finite(age),
      !is.na(sex),
      !is.na(project),
      !is.na(stage_major),
      !is.na(stage_group)
    )
  
  expression_mean <- mean(
    d[[gene]],
    na.rm = TRUE
  )
  
  expression_sd <- sd(
    d[[gene]],
    na.rm = TRUE
  )
  
  if (
    !is.finite(expression_sd) ||
    expression_sd <= 0
  ) {
    return(NULL)
  }
  
  d <- d |>
    mutate(
      expr_z =
        (.data[[gene]] - expression_mean) /
        expression_sd
    )
  
  reduced_model <- survival::coxph(
    survival::Surv(
      OS_days,
      event
    ) ~
      expr_z +
      age +
      sex +
      project +
      stage_major,
    data = d,
    x = TRUE,
    model = TRUE
  )
  
  full_model <- survival::coxph(
    survival::Surv(
      OS_days,
      event
    ) ~
      expr_z +
      expr_z:stage_group +
      age +
      sex +
      project +
      stage_major,
    data = d,
    x = TRUE,
    model = TRUE
  )
  
  sm <- summary(
    full_model
  )
  
  interaction_row <- grep(
    "expr_z:stage_groupStage_III_IV|stage_groupStage_III_IV:expr_z",
    rownames(
      sm$coefficients
    ),
    value = TRUE
  )
  
  if (length(interaction_row) != 1) {
    stop(
      "Interaction coefficient not found for ",
      gene
    )
  }
  
  interaction_row <- interaction_row[1]
  
  lrt <- anova(
    reduced_model,
    full_model,
    test = "LRT"
  )
  
  p_column <- grep(
    "Pr\\(",
    colnames(lrt),
    value = TRUE
  )
  
  if (length(p_column) == 0) {
    stop(
      "Likelihood-ratio p-value was not returned for ",
      gene
    )
  }
  
  interaction_lrt_p <- as.numeric(
    lrt[
      nrow(lrt),
      p_column[1]
    ]
  )
  
  covariance <- vcov(
    full_model
  )
  
  beta_early <- coef(
    full_model
  )[
    "expr_z"
  ]
  
  se_early <- sqrt(
    covariance[
      "expr_z",
      "expr_z"
    ]
  )
  
  beta_interaction <- coef(
    full_model
  )[
    interaction_row
  ]
  
  beta_advanced <-
    beta_early +
    beta_interaction
  
  variance_advanced <-
    covariance[
      "expr_z",
      "expr_z"
    ] +
    covariance[
      interaction_row,
      interaction_row
    ] +
    2 *
    covariance[
      "expr_z",
      interaction_row
    ]
  
  se_advanced <- sqrt(
    variance_advanced
  )
  
  ph <- tryCatch(
    survival::cox.zph(
      full_model,
      terms = TRUE
    ),
    error = function(e) NULL
  )
  
  ph_global <- NA_real_
  
  if (
    !is.null(ph) &&
    "GLOBAL" %in%
    rownames(ph$table)
  ) {
    ph_global <- ph$table[
      "GLOBAL",
      "p"
    ]
  }
  
  list(
    fit = full_model,
    summary = tibble(
      gene = gene,
      n = nrow(d),
      events = sum(
        d$event,
        na.rm = TRUE
      ),
      n_Stage_I_II = sum(
        d$stage_group ==
          "Stage_I_II"
      ),
      events_Stage_I_II = sum(
        d$event[
          d$stage_group ==
            "Stage_I_II"
        ],
        na.rm = TRUE
      ),
      n_Stage_III_IV = sum(
        d$stage_group ==
          "Stage_III_IV"
      ),
      events_Stage_III_IV = sum(
        d$event[
          d$stage_group ==
            "Stage_III_IV"
        ],
        na.rm = TRUE
      ),
      HR_Stage_I_II = exp(
        beta_early
      ),
      CI_low_Stage_I_II = exp(
        beta_early -
          1.96 * se_early
      ),
      CI_high_Stage_I_II = exp(
        beta_early +
          1.96 * se_early
      ),
      p_Stage_I_II =
        2 *
        pnorm(
          abs(
            beta_early /
              se_early
          ),
          lower.tail = FALSE
        ),
      HR_Stage_III_IV = exp(
        beta_advanced
      ),
      CI_low_Stage_III_IV = exp(
        beta_advanced -
          1.96 * se_advanced
      ),
      CI_high_Stage_III_IV = exp(
        beta_advanced +
          1.96 * se_advanced
      ),
      p_Stage_III_IV =
        2 *
        pnorm(
          abs(
            beta_advanced /
              se_advanced
          ),
          lower.tail = FALSE
        ),
      interaction_HR_ratio =
        sm$coefficients[
          interaction_row,
          "exp(coef)"
        ],
      interaction_LRT_p =
        interaction_lrt_p,
      interaction_Wald_p =
        sm$coefficients[
          interaction_row,
          "Pr(>|z|)"
        ],
      PH_global_p =
        ph_global
    )
  )
}

stage_interaction_fits <- lapply(
  FOCAL_GENES,
  fit_stage_interaction,
  data = survival_data_tpm
)

names(stage_interaction_fits) <- FOCAL_GENES

stage_interaction_fits <- stage_interaction_fits[
  !vapply(
    stage_interaction_fits,
    is.null,
    logical(1)
  )
]

stage_interaction_summary <- bind_rows(
  lapply(
    stage_interaction_fits,
    `[[`,
    "summary"
  )
) |>
  mutate(
    interaction_FDR = p.adjust(
      interaction_LRT_p,
      method = "BH"
    )
  )

write_tsv_safe(
  stage_interaction_summary,
  "TCGA_TPM_robust_stage_interaction.tsv"
)

# Publication-quality stage-interaction forest plot is rendered in Section 28.


###############################################################################
# 25. KAPLAN-MEIER SURVIVAL CURVES
###############################################################################

make_stage_km <- function(
    gene,
    group_value,
    group_label,
    panel_label,
    data
) {
  
  d <- data |>
    filter(
      stage_group == group_value,
      !is.na(.data[[gene]]),
      is.finite(.data[[gene]])
    )
  
  if (nrow(d) < 2) {
    return(NULL)
  }
  
  cutoff <- median(
    d[[gene]],
    na.rm = TRUE
  )
  
  d <- d |>
    mutate(
      expression_group = factor(
        if_else(
          .data[[gene]] >= cutoff,
          "High",
          "Low"
        ),
        levels = c(
          "Low",
          "High"
        )
      )
    )
  
  fit <- survival::survfit(
    survival::Surv(
      OS_days,
      event
    ) ~
      expression_group,
    data = d
  )
  
  logrank <- survival::survdiff(
    survival::Surv(
      OS_days,
      event
    ) ~
      expression_group,
    data = d
  )
  
  logrank_p <- pchisq(
    logrank$chisq,
    df = 1,
    lower.tail = FALSE
  )
  
  binary_cox <- survival::coxph(
    survival::Surv(
      OS_days,
      event
    ) ~
      expression_group,
    data = d
  )
  
  binary_summary <- summary(
    binary_cox
  )
  
  binary_ci <- exp(
    confint(
      binary_cox
    )[
      "expression_groupHigh",
    ]
  )
  
  summary_row <- tibble(
    gene = gene,
    stage_group = group_label,
    n = nrow(d),
    events = sum(
      d$event,
      na.rm = TRUE
    ),
    median_cutoff = cutoff,
    n_low = sum(
      d$expression_group ==
        "Low"
    ),
    n_high = sum(
      d$expression_group ==
        "High"
    ),
    events_low = sum(
      d$event[
        d$expression_group ==
          "Low"
      ],
      na.rm = TRUE
    ),
    events_high = sum(
      d$event[
        d$expression_group ==
          "High"
      ],
      na.rm = TRUE
    ),
    HR_high_vs_low =
      binary_summary$coefficients[
        "expression_groupHigh",
        "exp(coef)"
      ],
    CI_low = binary_ci[1],
    CI_high = binary_ci[2],
    logrank_p = logrank_p
  )
  
  plot_object <- survminer::ggsurvplot(
    fit,
    data = d,
    pval = TRUE,
    pval.method = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.24,
    risk.table.col = "strata",
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      unname(PLOT_COLORS["low"]),
      unname(PLOT_COLORS["high"])
    ),
    xlab = "Time (days)",
    ylab = "Overall survival probability",
    title = paste0(
      panel_label,
      "  ",
      gene,
      " | ",
      group_label
    ),
    legend.title = "Expression",
    legend.labs = c(
      "Low",
      "High"
    ),
    ggtheme = theme_classic(
      base_size = 11
    ),
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE
  )
  
  list(
    summary = summary_row,
    plot = plot_object
  )
}

km_definitions <- tribble(
  ~gene,    ~stage_group,    ~stage_label,     ~panel,
  "AQP8",   "Stage_I_II",   "Stage I-II",     "A",
  "GUCA2A", "Stage_I_II",   "Stage I-II",     "B",
  "MS4A12", "Stage_I_II",   "Stage I-II",     "C",
  "AQP8",   "Stage_III_IV", "Stage III-IV",   "D",
  "GUCA2A", "Stage_III_IV", "Stage III-IV",   "E",
  "MS4A12", "Stage_III_IV", "Stage III-IV",   "F"
)

km_results <- list()
km_plots <- list()

for (i in seq_len(nrow(km_definitions))) {
  
  row <- km_definitions[
    i,
  ]
  
  result <- make_stage_km(
    gene = row$gene,
    group_value = row$stage_group,
    group_label = row$stage_label,
    panel_label = row$panel,
    data = survival_data_tpm
  )
  
  if (is.null(result)) {
    next
  }
  
  key <- paste(
    row$stage_group,
    row$gene,
    sep = "__"
  )
  
  km_results[[key]] <- result$summary
  km_plots[[key]] <- result$plot
}

km_summary <- bind_rows(
  km_results
) |>
  group_by(
    stage_group
  ) |>
  mutate(
    FDR_three_genes = p.adjust(
      logrank_p,
      method = "BH"
    )
  ) |>
  ungroup()

write_tsv_safe(
  km_summary,
  "TCGA_TPM_KM_stage_median.tsv"
)

# Publication-quality six-panel Kaplan-Meier figure is rendered in Section 28.

###############################################################################
# 26. INTEGRATIVE OVERLAP SUMMARY
###############################################################################

# Symbols reused by the integrative overlap summary, key-count table, and
# publication overlap figures. Define them before their first use so a clean
# sequential run does not depend on objects created later in Section 28.
up_symbols <- unique(
  stats::na.omit(
    meta_up$symbol
  )
)

down_symbols <- unique(
  stats::na.omit(
    meta_down$symbol
  )
)

turquoise_hub_symbols <- unique(
  stats::na.omit(
    turquoise_hubs$symbol
  )
)

blue_hub_symbols <- unique(
  stats::na.omit(
    blue_hubs$symbol
  )
)

mcode4 <- mcode_modules[["4"]]

tcga_top25_down_symbols <- unique(
  na.omit(
    tcga_coad_read_top25_down$symbol
  )
)

candidate_pool <- unique(
  c(
    FOCAL_GENES,
    intersect(
      down_symbols,
      blue_hub_symbols
    ),
    mcode4
  )
)

ppi_candidate_annotation <- NULL

if (!is.null(ppi_centrality)) {
  
  ppi_annotation_cols <- c(
    "symbol",
    unlist(
      lapply(
        ppi_metric_cols,
        function(metric) {
          c(
            metric,
            paste0(
              metric,
              "_rank"
            ),
            paste0(
              metric,
              "_quartile"
            ),
            paste0(
              metric,
              "_Q1_or_Q2"
            )
          )
        }
      )
    ),
    "final_score",
    "final_score_rank",
    "final_score_quartile",
    "final_score_Q1_or_Q2",
    "n_metrics_Q1",
    "n_metrics_Q1_or_Q2"
  )
  
  ppi_candidate_annotation <- ppi_centrality |>
    dplyr::select(
      any_of(
        ppi_annotation_cols
      )
    )
}

candidate_summary <- tibble(
  symbol = candidate_pool
) |>
  mutate(
    meta_downregulated =
      symbol %in%
      down_symbols,
    
    normal_module_hub =
      symbol %in%
      blue_hub_symbols,
    
    MCODE_module4 =
      symbol %in%
      mcode4,
    
    TCGA_COAD_READ_top25_down =
      symbol %in%
      tcga_top25_down_symbols
  ) |>
  left_join(
    meta_results |>
      dplyr::select(
        symbol,
        meta_logFC,
        FDR,
        I2
      ) |>
      distinct(
        symbol,
        .keep_all = TRUE
      ),
    by = "symbol"
  )

if (!is.null(ppi_candidate_annotation)) {
  candidate_summary <- candidate_summary |>
    left_join(
      ppi_candidate_annotation,
      by = "symbol"
    )
}

# Additional PPI evidence used in Figure S12:
# TRUE means the gene lies in Q1 or Q2 of the network-wide Final Centrality
# ranking, i.e. above the median. Genes absent from the reconstructed STRING
# network are coded FALSE for this presence/absence evidence panel.
if ("final_score_Q1_or_Q2" %in% colnames(candidate_summary)) {
  candidate_summary <- candidate_summary |>
    mutate(
      final_centrality_above_median = dplyr::coalesce(
        as.logical(final_score_Q1_or_Q2),
        FALSE
      )
    )
} else {
  candidate_summary$final_centrality_above_median <- FALSE
}

candidate_summary <- candidate_summary |>
  left_join(
    overall_cox_summary |>
      dplyr::select(
        gene,
        HR_per_SD,
        CI_low,
        CI_high,
        p_value,
        FDR
      ) |>
      dplyr::rename(
        symbol = gene,
        Overall_Cox_HR = HR_per_SD,
        Overall_Cox_CI_low = CI_low,
        Overall_Cox_CI_high = CI_high,
        Overall_Cox_p = p_value,
        Overall_Cox_FDR = FDR
      ),
    by = "symbol"
  ) |>
  left_join(
    stage_interaction_summary |>
      dplyr::select(
        gene,
        HR_Stage_I_II,
        p_Stage_I_II,
        HR_Stage_III_IV,
        p_Stage_III_IV,
        interaction_LRT_p,
        interaction_FDR
      ) |>
      dplyr::rename(
        symbol = gene
      ),
    by = "symbol"
  ) |>
  arrange(
    desc(
      meta_downregulated
    ),
    desc(
      normal_module_hub
    ),
    desc(
      MCODE_module4
    ),
    FDR
  )

write_tsv_safe(
  candidate_summary,
  "Integrated_candidate_gene_summary.tsv"
)

evidence_columns <- c(
  "meta_downregulated",
  "normal_module_hub",
  "MCODE_module4",
  "final_centrality_above_median",
  "TCGA_COAD_READ_top25_down"
)

integrative_plot_df <- candidate_summary |>
  dplyr::select(
    symbol,
    all_of(
      evidence_columns
    )
  ) |>
  mutate(
    evidence_count = rowSums(
      across(
        all_of(
          evidence_columns
        )
      ),
      na.rm = TRUE
    )
  ) |>
  arrange(
    evidence_count,
    symbol
  ) |>
  dplyr::select(
    -evidence_count
  ) |>
  pivot_longer(
    cols = all_of(
      evidence_columns
    ),
    names_to = "evidence",
    values_to = "present"
  ) |>
  mutate(
    symbol = factor(
      symbol,
      levels = unique(
        symbol
      )
    ),
    evidence = factor(
      evidence,
      levels = evidence_columns,
      labels = c(
        "Meta-analysis downregulated",
        "Normal-associated WGCNA hub",
        "MCODE module 4",
        "Final Centrality above median",
        "TCGA top-25 downregulated"
      )
    )
  )

p_integrative_evidence <- ggplot(
  integrative_plot_df,
  aes(
    x = evidence,
    y = symbol,
    fill = present
  )
) +
  geom_tile(
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      `FALSE` = "grey92",
      `TRUE` = unname(PLOT_COLORS["accent"])
    )
  ) +
  labs(
    title = "Integrated evidence across analysis layers",
    x = NULL,
    y = NULL,
    fill = "Present"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    axis.text.y = element_text(
      size = 7
    )
  )

save_supplementary_plot(
  p_integrative_evidence,
  "Figure_S12_Integrated_candidate_evidence",
  width = 10.5,
  height = min(
    40,
    max(
      6,
      0.24 *
        length(
          unique(
            candidate_summary$symbol
          )
        )
    )
  )
)

###############################################################################
# 27. KEY COUNTS
###############################################################################

n_geo_pairs <- sum(
  cohort_summary$n_pairs
)

overlap_summary <- tibble(
  metric = c(
    "Number of paired GEO patients",
    "Number of GEO samples",
    "Meta-analysis upregulated DEGs",
    "Meta-analysis downregulated DEGs",
    paste0(
      "Up DEGs in ",
      turquoise_module,
      " module"
    ),
    paste0(
      "Down DEGs in ",
      blue_module,
      " module"
    ),
    paste0(
      "Hub genes in ",
      turquoise_module,
      " module"
    ),
    paste0(
      "Hub genes in ",
      blue_module,
      " module"
    )
  ),
  value = c(
    n_geo_pairs,
    2 * n_geo_pairs,
    nrow(meta_up),
    nrow(meta_down),
    length(
      intersect(
        up_symbols,
        module_map |>
          filter(
            module ==
              turquoise_module
          ) |>
          pull(symbol)
      )
    ),
    length(
      intersect(
        down_symbols,
        module_map |>
          filter(
            module ==
              blue_module
          ) |>
          pull(symbol)
      )
    ),
    nrow(turquoise_hubs),
    nrow(blue_hubs)
  )
)

write_tsv_safe(
  overlap_summary,
  "Manuscript_key_counts_corrected.tsv"
)

###############################################################################
# 28. PUBLICATION FIGURES
###############################################################################

# Final rendering is centralized here so that analytical calculations and
# publication styling remain separate and reproducible.

publication_required_objects <- c(
  "fig1",
  "go_plot_df",
  "kegg_plot_df",
  "g",
  "mcode_modules",
  "ppi_centrality",
  "sft_df",
  "WGCNA_SOFT_POWER_USED",
  "module_assoc",
  "module_map",
  "turquoise_module",
  "blue_module",
  "turquoise_hubs",
  "blue_hubs",
  "tpm_all",
  "tpm_rows",
  "clinical_stage_data",
  "tpm_expression_patient",
  "stage_interaction_summary",
  "km_plots",
  "km_summary"
)

missing_publication_objects <- publication_required_objects[
  !vapply(
    publication_required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_publication_objects) > 0) {
  stop(
    "Publication figure rendering cannot continue. Missing object(s): ",
    paste(
      missing_publication_objects,
      collapse = ", "
    )
  )
}

if (
  is.null(
    ppi_centrality
  ) ||
  nrow(
    ppi_centrality
  ) == 0
) {
  stop(
    "Publication figure rendering requires a non-empty PPI centrality table."
  )
}


###############################################################################
# Figure 1 — meta-analysis differential expression
###############################################################################

p_volcano_pub <-
  p_volcano +
  ggplot2::labs(
    title = "Meta-analysis volcano plot"
  )

p_ma_pub <-
  p_ma +
  ggplot2::labs(
    title = "Meta-analysis MA plot"
  )

figure1_pub <-
  p_volcano_pub +
  p_ma_pub +
  patchwork::plot_layout(
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      plot.tag.position = "topleft"
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )


###############################################################################
# Figures 2 and 3 — enrichment styling
###############################################################################

down_low  <- "#DCEAF7"
down_high <- "#2166AC"

up_low    <- "#F7D9D5"
up_high   <- "#B2182B"

point_outline <- "grey25"


prepare_enrichment_df <- function(
    x
) {
  
  x |>
    dplyr::mutate(
      direction = factor(
        direction,
        levels = c(
          "Down",
          "Up"
        )
      ),
      neglog10_fdr = -log10(
        p.adjust
      )
    ) |>
    dplyr::filter(
      is.finite(
        GeneRatioNum
      ),
      is.finite(
        neglog10_fdr
      ),
      is.finite(
        Count
      )
    )
}


go_pub_df <- prepare_enrichment_df(
  go_plot_df
)

kegg_pub_df <- prepare_enrichment_df(
  kegg_plot_df
)


###############################################################################
# Helpers
###############################################################################

direction_fill_scale <- function(
    direction,
    legend_title
) {
  
  if (identical(direction, "Down")) {
    
    ggplot2::scale_fill_gradient(
      low = down_low,
      high = down_high,
      name = legend_title,
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(
          26,
          "mm"
        ),
        barheight = grid::unit(
          3.4,
          "mm"
        )
      )
    )
    
  } else {
    
    ggplot2::scale_fill_gradient(
      low = up_low,
      high = up_high,
      name = legend_title,
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = grid::unit(
          26,
          "mm"
        ),
        barheight = grid::unit(
          3.4,
          "mm"
        )
      )
    )
  }
}


common_size_scale <- function(
    limits_value
) {
  
  ggplot2::scale_size_continuous(
    name = "Genes",
    range = c(
      2.5,
      7.3
    ),
    limits = limits_value,
    breaks = scales::pretty_breaks(
      n = 3
    )
  )
}


panel_theme <- function(
    y_size
) {
  
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 5
        )
      ),
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      
      plot.tag.position = c(
        0.01,
        1.015
      ),
      
      axis.title.x = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(
          t = 7
        )
      ),
      
      axis.text.y = ggplot2::element_text(
        size = y_size,
        color = "grey15"
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 8.5,
        color = "grey20"
      ),
      
      strip.background = ggplot2::element_blank(),
      
      strip.text.y.right = ggplot2::element_text(
        face = "bold",
        angle = 0,
        size = 9.2,
        color = "grey20",
        margin = ggplot2::margin(
          l = 6
        )
      ),
      
      legend.position = "bottom",
      
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 8.5
      ),
      
      legend.text = ggplot2::element_text(
        size = 7.8
      ),
      
      legend.box = "horizontal",
      
      plot.margin = ggplot2::margin(
        t = 10,
        r = 10,
        b = 7,
        l = 8
      )
    )
}


###############################################################################
# Figure 2 — GO enrichment
###############################################################################

go_size_limits <- range(
  go_pub_df$Count,
  na.rm = TRUE
)

go_x_limits <- range(
  go_pub_df$GeneRatioNum,
  na.rm = TRUE
)

go_x_pad <- diff(
  go_x_limits
) * 0.06

if (
  !is.finite(go_x_pad) ||
  go_x_pad <= 0
) {
  go_x_pad <- max(
    go_x_limits,
    na.rm = TRUE
  ) * 0.05
}

go_x_display <- c(
  max(
    0,
    go_x_limits[1] - go_x_pad
  ),
  go_x_limits[2] + go_x_pad
)


make_go_direction_panel <- function(
    direction,
    tag,
    title_text
) {
  
  d <- go_pub_df |>
    dplyr::filter(
      direction == !!direction
    ) |>
    dplyr::mutate(
      ontology = factor(
        ontology,
        levels = c(
          "BP",
          "CC",
          "MF"
        )
      )
    )
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = GeneRatioNum,
      y = stats::reorder(
        Description,
        GeneRatioNum
      ),
      size = Count,
      fill = neglog10_fdr
    )
  ) +
    
    ggplot2::geom_point(
      shape = 21,
      color = point_outline,
      stroke = 0.24,
      alpha = 0.97
    ) +
    
    ggplot2::facet_grid(
      ontology ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    
    direction_fill_scale(
      direction = direction,
      legend_title = expression(
        -log[10](FDR)
      )
    ) +
    
    common_size_scale(
      go_size_limits
    ) +
    
    ggplot2::coord_cartesian(
      xlim = go_x_display,
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = title_text,
      tag = tag,
      x = "Gene ratio",
      y = NULL
    ) +
    
    panel_theme(
      y_size = 7.3
    )
}


go_down_panel <- make_go_direction_panel(
  direction = "Down",
  tag = "A",
  title_text = "Downregulated DEGs"
)

go_up_panel <- make_go_direction_panel(
  direction = "Up",
  tag = "B",
  title_text = "Upregulated DEGs"
)


figure2_pub <-
  go_down_panel +
  go_up_panel +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    ),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "GO enrichment of meta-analysis DEGs",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 8
        )
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )


###############################################################################
# Figure 3 — KEGG enrichment
###############################################################################

kegg_size_limits <- range(
  kegg_pub_df$Count,
  na.rm = TRUE
)

kegg_x_limits <- range(
  kegg_pub_df$GeneRatioNum,
  na.rm = TRUE
)

kegg_x_pad <- diff(
  kegg_x_limits
) * 0.06

if (
  !is.finite(kegg_x_pad) ||
  kegg_x_pad <= 0
) {
  kegg_x_pad <- max(
    kegg_x_limits,
    na.rm = TRUE
  ) * 0.05
}

kegg_x_display <- c(
  max(
    0,
    kegg_x_limits[1] - kegg_x_pad
  ),
  kegg_x_limits[2] + kegg_x_pad
)


make_kegg_direction_panel <- function(
    direction,
    tag,
    title_text
) {
  
  d <- kegg_pub_df |>
    dplyr::filter(
      direction == !!direction
    )
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = GeneRatioNum,
      y = stats::reorder(
        Description,
        GeneRatioNum
      ),
      size = Count,
      fill = neglog10_fdr
    )
  ) +
    
    ggplot2::geom_point(
      shape = 21,
      color = point_outline,
      stroke = 0.24,
      alpha = 0.97
    ) +
    
    direction_fill_scale(
      direction = direction,
      legend_title = expression(
        -log[10](FDR)
      )
    ) +
    
    common_size_scale(
      kegg_size_limits
    ) +
    
    ggplot2::coord_cartesian(
      xlim = kegg_x_display,
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = title_text,
      tag = tag,
      x = "Gene ratio",
      y = NULL
    ) +
    
    panel_theme(
      y_size = 8.0
    )
}


kegg_down_panel <- make_kegg_direction_panel(
  direction = "Down",
  tag = "A",
  title_text = "Downregulated DEGs"
)

kegg_up_panel <- make_kegg_direction_panel(
  direction = "Up",
  tag = "B",
  title_text = "Upregulated DEGs"
)


figure3_pub <-
  kegg_down_panel +
  kegg_up_panel +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    ),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "KEGG pathway enrichment of meta-analysis DEGs",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 8
        )
      )
    )
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )



###############################################################################
# Figure 4 — PPI/MCODE network styling
###############################################################################

effect_map <- stats::setNames(
  meta_results$meta_logFC,
  meta_results$symbol
)

score_map <- stats::setNames(
  ppi_centrality$final_score,
  ppi_centrality$symbol
)

all_module_genes <- unique(
  unlist(
    mcode_modules,
    use.names = FALSE
  )
)

all_module_genes <- intersect(
  all_module_genes,
  igraph::V(g)$name
)

centrality_values <- unname(
  score_map[
    all_module_genes
  ]
)

centrality_values <- centrality_values[
  is.finite(
    centrality_values
  )
]

centrality_limits <- range(
  centrality_values,
  na.rm = TRUE
)

if (
  length(centrality_limits) != 2 ||
  any(
    !is.finite(
      centrality_limits
    )
  ) ||
  diff(
    centrality_limits
  ) == 0
) {
  centrality_limits <- c(
    0,
    1
  )
}


###############################################################################
# Module network panels
###############################################################################

make_mcode_panel <- function(
    module_id,
    panel_tag
) {
  
  genes <- intersect(
    mcode_modules[[module_id]],
    igraph::V(g)$name
  )
  
  if (length(genes) < 2) {
    return(NULL)
  }
  
  sg <- igraph::induced_subgraph(
    g,
    vids = genes
  )
  
  igraph::V(sg)$effect <- unname(
    effect_map[
      igraph::V(sg)$name
    ]
  )
  
  igraph::V(sg)$centrality <- unname(
    score_map[
      igraph::V(sg)$name
    ]
  )
  
  # A fixed seed keeps the network layout reproducible.
  set.seed(
    4100 +
      as.integer(
        module_id
      )
  )
  
  label_size <- dplyr::case_when(
    igraph::vcount(sg) >= 30 ~ 2.15,
    igraph::vcount(sg) >= 15 ~ 2.40,
    TRUE ~ 2.65
  )
  
  edge_alpha <- 0.25
  
  ggraph::ggraph(
    sg,
    layout = "fr"
  ) +
    
    ggraph::geom_edge_link(
      color = "grey55",
      alpha = edge_alpha,
      linewidth = 0.40,
      show.legend = FALSE
    ) +
    
    ggraph::geom_node_point(
      ggplot2::aes(
        fill = effect,
        size = centrality
      ),
      shape = 21,
      color = "grey20",
      stroke = 0.30
    ) +
    
    ggraph::geom_node_text(
      ggplot2::aes(
        label = name
      ),
      repel = TRUE,
      size = label_size,
      color = "grey12",
      box.padding = 0.42,
      point.padding = 0.24,
      force = 1.25,
      max.overlaps = Inf,
      max.iter = 10000,
      min.segment.length = 0,
      segment.color = "grey65",
      segment.alpha = 0.70,
      segment.size = 0.25,
      show.legend = FALSE
    ) +
    
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    
    ggplot2::scale_size_continuous(
      range = c(
        2.7,
        6.5
      ),
      limits = centrality_limits
    ) +
    
    ggplot2::labs(
      title = paste0(
        "MCODE module ",
        module_id
      ),
      tag = panel_tag,
      fill = "Meta log2FC",
      size = "Final Centrality score"
    ) +
    
    ggplot2::theme_void(
      base_size = 11
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 5
        )
      ),
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      
      plot.tag.position = c(
        0.015,
        0.985
      ),
      
      plot.margin = ggplot2::margin(
        t = 12,
        r = 18,
        b = 12,
        l = 18
      ),
      
      legend.position = "none"
    )
}


module_ids <- names(
  mcode_modules
)

module_ids <- module_ids[
  module_ids %in%
    as.character(
      seq_len(
        min(
          5,
          length(
            module_ids
          )
        )
      )
    )
]

panel_tags <- LETTERS[
  seq_along(
    module_ids
  )
]

module_plots_pub <- stats::setNames(
  lapply(
    seq_along(
      module_ids
    ),
    function(i) {
      make_mcode_panel(
        module_id = module_ids[i],
        panel_tag = panel_tags[i]
      )
    }
  ),
  module_ids
)

module_plots_pub <- module_plots_pub[
  !vapply(
    module_plots_pub,
    is.null,
    logical(1)
  )
]

if (length(module_plots_pub) < 5) {
  stop(
    "Five MCODE module panels are required for Figure 4."
  )
}


# A compact custom legend is shown only once, below module 5.
# Node colors are scaled within each module, so the color key is intentionally
# qualitative and shows the full blue-white-red direction of Meta log2FC.

legend_gradient <- data.frame(
  x = seq(
    0.65,
    4.35,
    length.out = 121
  ),
  value = seq(
    -1,
    1,
    length.out = 121
  )
)

figure4_legend <-
  ggplot2::ggplot() +
  
  ggplot2::geom_tile(
    data = legend_gradient,
    ggplot2::aes(
      x = x,
      y = 1.05,
      fill = value
    ),
    width = 0.035,
    height = 0.34
  ) +
  
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    guide = "none"
  ) +
  
  ggplot2::annotate(
    "text",
    x = 2.50,
    y = 1.78,
    label = "Meta log2FC",
    fontface = "bold",
    size = 3.15
  ) +
  
  ggplot2::annotate(
    "text",
    x = c(
      0.65,
      2.50,
      4.35
    ),
    y = 0.53,
    label = c(
      "Downregulated",
      "0",
      "Upregulated"
    ),
    size = 2.55,
    color = "grey25"
  ) +
  
  ggplot2::annotate(
    "text",
    x = 7.35,
    y = 1.78,
    label = "Final Centrality score",
    fontface = "bold",
    size = 3.15
  ) +
  
  ggplot2::geom_point(
    data = data.frame(
      x = c(
        6.35,
        7.35,
        8.35
      ),
      y = rep(
        1.05,
        3
      ),
      point_size = c(
        3.0,
        4.8,
        6.5
      )
    ),
    ggplot2::aes(
      x = x,
      y = y,
      size = point_size
    ),
    shape = 21,
    fill = "grey85",
    color = "grey25",
    stroke = 0.35,
    show.legend = FALSE
  ) +
  
  ggplot2::scale_size_identity() +
  
  ggplot2::annotate(
    "text",
    x = c(
      6.35,
      8.35
    ),
    y = 0.53,
    label = c(
      "Lower",
      "Higher"
    ),
    size = 2.55,
    color = "grey25"
  ) +
  
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      9
    ),
    ylim = c(
      0.25,
      2
    ),
    clip = "off"
  ) +
  
  ggplot2::theme_void() +
  
  ggplot2::theme(
    plot.margin = ggplot2::margin(
      t = 0,
      r = 4,
      b = 0,
      l = 4
    )
  )


###############################################################################
# Figure assembly
###############################################################################

row_1 <-
  module_plots_pub[["1"]] +
  module_plots_pub[["2"]] +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )


row_2 <-
  module_plots_pub[["3"]] +
  module_plots_pub[["4"]] +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )


# The fifth module is centered deliberately so the final row does not look
# unfinished or visually heavier on the left.
module_5_with_legend <-
  module_plots_pub[["5"]] /
  figure4_legend +
  patchwork::plot_layout(
    heights = c(
      1,
      0.24
    )
  )


row_3 <-
  patchwork::plot_spacer() +
  module_5_with_legend +
  patchwork::plot_spacer() +
  patchwork::plot_layout(
    widths = c(
      0.50,
      1,
      0.50
    )
  )


figure4_pub <-
  row_1 /
  row_2 /
  row_3 +
  patchwork::plot_layout(
    heights = c(
      1,
      1,
      1.24
    )
  )



# Figure 5 palette
col_blue   <- unname(PLOT_COLORS["down"])
col_red    <- unname(PLOT_COLORS["up"])
col_accent <- unname(PLOT_COLORS["accent"])
col_grey   <- unname(PLOT_COLORS["neutral"])

###############################################################################
# Figure 5 — WGCNA diagnostics, gene dendrogram, and module association
###############################################################################

selected_sft <- sft_df |>
  dplyr::filter(
    Power == WGCNA_SOFT_POWER_USED
  )

if (nrow(selected_sft) != 1) {
  stop("The selected soft-threshold power could not be identified uniquely.")
}


# Panel A: scale-free topology fit

p_sft1_pub <- ggplot2::ggplot(
  sft_df,
  ggplot2::aes(
    x = Power,
    y = SFT.R.sq
  )
) +
  ggplot2::geom_line(
    color = scales::alpha(
      col_blue,
      0.62
    ),
    linewidth = 0.75
  ) +
  ggplot2::geom_point(
    color = col_blue,
    fill = "white",
    shape = 21,
    stroke = 0.8,
    size = 2.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = Power
    ),
    color = "grey35",
    size = 2.75,
    nudge_y = 0.035,
    check_overlap = FALSE
  ) +
  ggplot2::geom_hline(
    yintercept = WGCNA_SFT_R2_TARGET,
    color = "grey45",
    linetype = "dashed",
    linewidth = 0.55
  ) +
  ggplot2::geom_vline(
    xintercept = WGCNA_SOFT_POWER_USED,
    color = col_red,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    data = selected_sft,
    color = col_red,
    fill = col_red,
    shape = 21,
    stroke = 0.9,
    size = 4.2
  ) +
  ggplot2::geom_label(
    data = selected_sft,
    ggplot2::aes(
      label = paste0(
        "Selected power = ",
        Power
      )
    ),
    color = "white",
    fill = col_red,
    fontface = "bold",
    size = 3.05,
    label.size = 0,
    label.padding = grid::unit(
      0.16,
      "lines"
    ),
    nudge_x = 3.2,
    nudge_y = -0.095,
    show.legend = FALSE
  ) +
  ggplot2::annotate(
    "text",
    x = min(
      sft_df$Power,
      na.rm = TRUE
    ),
    y = WGCNA_SFT_R2_TARGET,
    label = paste0(
      " Target ",
      format(
        WGCNA_SFT_R2_TARGET,
        trim = TRUE
      )
    ),
    color = "grey35",
    hjust = 0,
    vjust = -0.55,
    size = 2.9
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(
      0,
      max(
        sft_df$Power,
        na.rm = TRUE
      ),
      by = 10
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.025,
        0.08
      )
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      max(
        1,
        sft_df$SFT.R.sq,
        na.rm = TRUE
      )
    ),
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.06
      )
    )
  ) +
  ggplot2::labs(
    title = "Scale-free topology fit",
    x = "Soft-threshold power",
    y = expression("Signed " * R^2)
  ) +
  ggplot2::theme_classic(
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = "topleft"
  ) +
  ggplot2::labs(tag = "A")


# Panel B: mean connectivity

connectivity_max <- max(
  sft_df$mean.k.,
  na.rm = TRUE
)

connectivity_bottom_pad <- 0.055 * connectivity_max
connectivity_top_pad <- 0.08 * connectivity_max

p_sft2_pub <- ggplot2::ggplot(
  sft_df,
  ggplot2::aes(
    x = Power,
    y = mean.k.
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "grey70",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_line(
    color = scales::alpha(
      col_blue,
      0.62
    ),
    linewidth = 0.75
  ) +
  ggplot2::geom_point(
    color = col_blue,
    fill = "white",
    shape = 21,
    stroke = 0.8,
    size = 2.5
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = Power
    ),
    color = "grey35",
    size = 2.65,
    nudge_y = 0.030 * connectivity_max,
    check_overlap = FALSE
  ) +
  ggplot2::geom_vline(
    xintercept = WGCNA_SOFT_POWER_USED,
    color = col_red,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::geom_point(
    data = selected_sft,
    color = col_red,
    fill = col_red,
    shape = 21,
    stroke = 0.9,
    size = 4.2
  ) +
  ggplot2::geom_label(
    data = selected_sft,
    ggplot2::aes(
      label = paste0(
        "Selected power = ",
        Power
      )
    ),
    color = "white",
    fill = col_red,
    fontface = "bold",
    size = 3.05,
    label.size = 0,
    label.padding = grid::unit(
      0.16,
      "lines"
    ),
    nudge_x = 3.0,
    nudge_y = 0.095 * connectivity_max,
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(
      0,
      max(
        sft_df$Power,
        na.rm = TRUE
      ),
      by = 10
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.025,
        0.08
      )
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = scales::breaks_pretty(
      n = 4
    ),
    expand = c(
      0,
      0
    )
  ) +
  ggplot2::coord_cartesian(
    ylim = c(
      -connectivity_bottom_pad,
      connectivity_max + connectivity_top_pad
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Mean connectivity",
    x = "Soft-threshold power",
    y = "Mean connectivity"
  ) +
  ggplot2::theme_classic(
    base_size = 11
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    axis.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = "topleft"
  ) +
  ggplot2::labs(tag = "B")


###############################################################################
# Module-trait association
###############################################################################

module_assoc_plot <- module_assoc |>
  dplyr::arrange(partial_r) |>
  dplyr::mutate(
    module = factor(
      module,
      levels = unique(module)
    ),
    module_label = as.character(module),
    module_text_color = dplyr::case_when(
      module_label == "yellow" ~ "#B8860B",
      module_label == "grey" ~ "#666666",
      module_label == "turquoise" ~ "#008C95",
      TRUE ~ module_label
    )
  )


p_modtrait_pub <- ggplot2::ggplot(
  module_assoc_plot,
  ggplot2::aes(
    x = 1,
    y = module,
    fill = partial_r
  )
) +
  ggplot2::geom_tile(
    width = 0.78,
    height = 0.95,
    color = "white",
    linewidth = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf(
        "r = %.2f\nFDR = %.2g",
        partial_r,
        FDR
      )
    ),
    size = 3.15,
    lineheight = 0.94
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 0.54,
      label = module_label,
      color = module_text_color
    ),
    hjust = 1,
    fontface = "bold",
    size = 3.7,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient2(
    low = col_blue,
    mid = "white",
    high = col_red,
    midpoint = 0
  ) +
  ggplot2::scale_x_continuous(
    breaks = NULL,
    expand = ggplot2::expansion(
      mult = 0
    )
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      0.38,
      1.48
    ),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Paired module–tissue association",
    x = NULL,
    y = NULL,
    fill = "Partial r"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(
        38,
        "mm"
      ),
      barheight = grid::unit(
        3.8,
        "mm"
      )
    )
  ) +
  ggplot2::theme_classic(
    base_size = 11
  ) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.margin = ggplot2::margin(
      t = 7,
      r = 10,
      b = 4,
      l = 36
    ),
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 13
    ),
    plot.tag.position = "topleft"
  ) +
  ggplot2::labs(tag = "D")


###############################################################################
# Figure assembly
###############################################################################

# Panel C: convert the base-R WGCNA dendrogram to a raster-backed ggplot.
# This avoids the device/capture problem that can cause Panels C and D to
# disappear when a base-graphics formula is embedded directly in patchwork.
dendro_tmp_png <- tempfile(
  pattern = "figure5_dendrogram_",
  fileext = ".png"
)

grDevices::png(
  filename = dendro_tmp_png,
  width = 5200,
  height = 2300,
  res = 400,
  bg = "white"
)

draw_wgcna_dendrogram()
grDevices::dev.off()

dendro_raster <- png::readPNG(
  dendro_tmp_png
)

unlink(
  dendro_tmp_png
)

p_dendro_pub <- ggplot2::ggplot() +
  ggplot2::annotation_custom(
    grob = grid::rasterGrob(
      dendro_raster,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc"),
      interpolate = TRUE
    ),
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE,
    clip = "off"
  ) +
  ggplot2::theme_void() +
  ggplot2::labs(
    tag = "C"
  ) +
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      face = "bold",
      size = 13,
      color = "black"
    ),
    plot.tag.position = "topleft",
    plot.margin = ggplot2::margin(
      t = 3,
      r = 3,
      b = 3,
      l = 3
    )
  )


# Use an explicit four-area design instead of nested patchworks.
# This guarantees that all four objects occupy their own defined region.
figure5_design <- c(
  patchwork::area(
    t = 1,
    l = 1,
    b = 1,
    r = 6
  ),
  patchwork::area(
    t = 1,
    l = 7,
    b = 1,
    r = 12
  ),
  patchwork::area(
    t = 2,
    l = 1,
    b = 2,
    r = 12
  ),
  patchwork::area(
    t = 3,
    l = 3,
    b = 3,
    r = 10
  )
)


# Final four-panel Figure 5:
# A = scale-free topology fit
# B = mean connectivity
# C = gene dendrogram and module colors
# D = paired module–tissue association
figure5_wgcna_pub <- patchwork::wrap_plots(
  p_sft1_pub,
  p_sft2_pub,
  p_dendro_pub,
  p_modtrait_pub,
  design = figure5_design
) +
  patchwork::plot_layout(
    heights = c(
      1.00,
      1.10,
      1.10
    )
  )

# Figure 6 palette
hub_colour <- "#B2182B"
nonhub_colour <- "grey72"
up_colour <- "#7A001F"
down_colour <- "#163A5F"

###############################################################################
# Figure 6 — WGCNA hubs and DEG overlap
###############################################################################

hub_scatter_pub <- function(
    module_name,
    label_genes = FOCAL_GENES
) {
  
  df <- module_map |>
    dplyr::filter(
      module == module_name,
      is.finite(MM),
      is.finite(GS)
    ) |>
    dplyr::mutate(
      abs_MM = abs(MM),
      label = dplyr::if_else(
        symbol %in% label_genes,
        symbol,
        NA_character_
      ),
      hub = factor(
        hub,
        levels = c(
          FALSE,
          TRUE
        )
      )
    )
  
  if (nrow(df) == 0) {
    stop(
      "No genes were available for module ",
      module_name,
      "."
    )
  }
  
  module_colour <- module_name
  
  if (module_name == "yellow") {
    module_colour <- "#B8860B"
  }
  
  if (module_name == "grey") {
    module_colour <- "#666666"
  }
  
  if (module_name == "turquoise") {
    module_colour <- "#008C95"
  }
  
  x_min <- min(
    df$abs_MM,
    na.rm = TRUE
  )
  
  x_max <- max(
    df$abs_MM,
    na.rm = TRUE
  )
  
  y_min <- min(
    df$GS,
    na.rm = TRUE
  )
  
  y_max <- max(
    df$GS,
    na.rm = TRUE
  )
  
  x_range <- x_max - x_min
  y_range <- y_max - y_min
  
  if (
    !is.finite(x_range) ||
    x_range <= 0
  ) {
    x_range <- 1
  }
  
  if (
    !is.finite(y_range) ||
    y_range <= 0
  ) {
    y_range <- 1
  }
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = abs_MM,
      y = GS
    )
  ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        color = hub
      ),
      alpha = 0.72,
      size = 1.55
    ) +
    
    ggplot2::scale_color_manual(
      values = c(
        `FALSE` = nonhub_colour,
        `TRUE` = hub_colour
      ),
      labels = c(
        `FALSE` = "FALSE",
        `TRUE` = "TRUE"
      ),
      drop = FALSE
    ) +
    
    ggplot2::geom_vline(
      xintercept = WGCNA_MM_CUTOFF,
      color = "grey45",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    
    ggplot2::geom_hline(
      yintercept = WGCNA_GS_CUTOFF,
      color = "grey45",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    
    ggplot2::annotate(
      "label",
      x = WGCNA_MM_CUTOFF,
      y = y_min + 0.035 * y_range,
      label = paste0(
        "|MM| = ",
        format(
          WGCNA_MM_CUTOFF,
          trim = TRUE
        )
      ),
      angle = 90,
      hjust = 0,
      vjust = -0.25,
      size = 2.8,
      label.size = 0,
      fill = scales::alpha(
        "white",
        0.88
      ),
      color = "grey30"
    ) +
    
    ggplot2::annotate(
      "label",
      x = x_min + 0.025 * x_range,
      y = WGCNA_GS_CUTOFF,
      label = paste0(
        "GS = ",
        format(
          WGCNA_GS_CUTOFF,
          trim = TRUE
        )
      ),
      hjust = 0,
      vjust = -0.35,
      size = 2.8,
      label.size = 0,
      fill = scales::alpha(
        "white",
        0.88
      ),
      color = "grey30"
    ) +
    
    ggrepel::geom_label_repel(
      data = df |>
        dplyr::filter(
          !is.na(label)
        ),
      ggplot2::aes(
        label = label
      ),
      color = hub_colour,
      fill = "white",
      fontface = "bold",
      size = 3.05,
      label.size = 0.25,
      box.padding = 0.35,
      point.padding = 0.20,
      min.segment.length = 0,
      segment.color = hub_colour,
      segment.alpha = 0.65,
      max.overlaps = Inf,
      seed = 1234,
      show.legend = FALSE
    ) +
    
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.04,
          0.08
        )
      )
    ) +
    
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0.05,
          0.08
        )
      )
    ) +
    
    ggplot2::labs(
      title = paste0(
        "Module ",
        module_name
      ),
      subtitle = paste0(
        "Hub genes: |MM| \u2265 ",
        format(
          WGCNA_MM_CUTOFF,
          trim = TRUE
        ),
        " and GS \u2265 ",
        format(
          WGCNA_GS_CUTOFF,
          trim = TRUE
        )
      ),
      x = "|Module membership|",
      y = "Paired gene significance |partial r|",
      color = "Hub"
    ) +
    
    ggplot2::coord_cartesian(
      clip = "off"
    ) +
    
    ggplot2::theme_classic(
      base_size = 11
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        color = module_colour
      ),
      plot.subtitle = ggplot2::element_text(
        size = 9.5,
        color = "grey35"
      ),
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "right",
      legend.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.key.height = grid::unit(
        4,
        "mm"
      ),
      plot.margin = ggplot2::margin(
        8,
        10,
        8,
        8
      )
    )
}


p_hub_t_pub <- hub_scatter_pub(
  turquoise_module
)

p_hub_b_pub <- hub_scatter_pub(
  blue_module
)


###############################################################################
# DEG overlap
###############################################################################

# up_symbols, down_symbols, turquoise_hub_symbols, and blue_hub_symbols were
# defined in Section 26 before their first use.

circle_points <- function(
    cx,
    cy,
    radius = 1,
    n = 400
) {
  
  angle <- seq(
    0,
    2 * pi,
    length.out = n
  )
  
  data.frame(
    x = cx + radius * cos(angle),
    y = cy + radius * sin(angle)
  )
}


make_two_set_venn <- function(
    set_top,
    set_bottom,
    label_top,
    label_bottom,
    title,
    colour_top,
    colour_bottom
) {
  
  set_top <- unique(
    stats::na.omit(
      set_top
    )
  )
  
  set_bottom <- unique(
    stats::na.omit(
      set_bottom
    )
  )
  
  only_top <- length(
    setdiff(
      set_top,
      set_bottom
    )
  )
  
  only_bottom <- length(
    setdiff(
      set_bottom,
      set_top
    )
  )
  
  overlap <- length(
    intersect(
      set_top,
      set_bottom
    )
  )
  
  union_n <- length(
    union(
      set_top,
      set_bottom
    )
  )
  
  pct <- function(n) {
    if (union_n == 0) {
      0
    } else {
      round(
        100 * n / union_n
      )
    }
  }
  
  # Slightly greater vertical separation makes the three numeric regions
  # visually distinct without changing any counts.
  top_circle <- circle_points(
    0,
    0.50,
    radius = 1.14
  )
  
  bottom_circle <- circle_points(
    0,
    -0.50,
    radius = 1.14
  )
  
  ggplot2::ggplot() +
    
    ggplot2::geom_polygon(
      data = top_circle,
      ggplot2::aes(
        x = x,
        y = y
      ),
      fill = scales::alpha(
        colour_top,
        0.68
      ),
      color = colour_top,
      linewidth = 0.85
    ) +
    
    ggplot2::geom_polygon(
      data = bottom_circle,
      ggplot2::aes(
        x = x,
        y = y
      ),
      fill = scales::alpha(
        colour_bottom,
        0.50
      ),
      color = colour_bottom,
      linewidth = 0.85
    ) +
    
    # Set names are outside the circle borders so they cannot overlap the arcs.
    ggplot2::annotate(
      "text",
      x = 0,
      y = 1.82,
      label = label_top,
      fontface = "bold",
      size = 3.8,
      color = colour_top
    ) +
    
    ggplot2::annotate(
      "text",
      x = 0,
      y = -1.82,
      label = label_bottom,
      fontface = "bold",
      size = 3.8,
      color = colour_bottom
    ) +
    
    # Counts are moved deeper into their own regions, away from circle borders.
    ggplot2::annotate(
      "text",
      x = 0,
      y = 1.02,
      label = sprintf(
        "%d
(%d%%)",
        only_top,
        pct(
          only_top
        )
      ),
      fontface = "bold",
      size = 3.75,
      lineheight = 0.96,
      color = "black"
    ) +
    
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = sprintf(
        "%d
(%d%%)",
        overlap,
        pct(
          overlap
        )
      ),
      fontface = "bold",
      size = 3.75,
      lineheight = 0.96,
      color = "black"
    ) +
    
    ggplot2::annotate(
      "text",
      x = 0,
      y = -1.02,
      label = sprintf(
        "%d
(%d%%)",
        only_bottom,
        pct(
          only_bottom
        )
      ),
      fontface = "bold",
      size = 3.75,
      lineheight = 0.96,
      color = "black"
    ) +
    
    ggplot2::coord_fixed(
      xlim = c(
        -1.42,
        1.42
      ),
      ylim = c(
        -2.02,
        2.02
      ),
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = title
    ) +
    
    ggplot2::theme_void(
      base_size = 11
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 4
        )
      ),
      plot.margin = ggplot2::margin(
        t = 8,
        r = 14,
        b = 8,
        l = 14
      )
    )
}


p_venn_up_pub <- make_two_set_venn(
  set_top = up_symbols,
  set_bottom = turquoise_hub_symbols,
  label_top = "Upregulated DEGs",
  label_bottom = "Turquoise hubs",
  title = "Upregulated overlap",
  colour_top = up_colour,
  colour_bottom = "turquoise"
)


p_venn_down_pub <- make_two_set_venn(
  set_top = down_symbols,
  set_bottom = blue_hub_symbols,
  label_top = "Downregulated DEGs",
  label_bottom = "Blue hubs",
  title = "Downregulated overlap",
  colour_top = down_colour,
  colour_bottom = "blue"
)


###############################################################################
# Assemble Figure 6
###############################################################################

top_row <- p_hub_t_pub + p_hub_b_pub +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    ),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "right"
  )


bottom_row <- p_venn_up_pub + p_venn_down_pub +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )


figure6_wgcna_pub <- top_row / bottom_row +
  patchwork::plot_layout(
    heights = c(
      1.10,
      0.90
    )
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      plot.tag.position = "topleft"
    )
  )



format_fdr <- function(x) {
  
  dplyr::case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ formatC(
      x,
      format = "e",
      digits = 2
    ),
    TRUE ~ formatC(
      x,
      format = "f",
      digits = 3
    )
  )
}


comparison_annotation <- function(
    data,
    group_var
) {
  
  data |>
    dplyr::group_by(
      gene
    ) |>
    dplyr::group_modify(
      ~{
        
        x <- .x
        group_values <- droplevels(
          factor(
            x[[group_var]]
          )
        )
        
        group_levels <- levels(
          group_values
        )
        
        if (length(group_levels) != 2L) {
          return(
            tibble::tibble(
              p_value = NA_real_
            )
          )
        }
        
        p_value <- tryCatch(
          stats::wilcox.test(
            x$expression[
              group_values == group_levels[1]
            ],
            x$expression[
              group_values == group_levels[2]
            ],
            exact = FALSE
          )$p.value,
          error = function(e) {
            NA_real_
          }
        )
        
        tibble::tibble(
          p_value = p_value
        )
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      FDR = stats::p.adjust(
        p_value,
        method = "BH"
      ),
      label = paste0(
        "FDR = ",
        vapply(
          FDR,
          format_fdr,
          character(1)
        )
      )
    )
}


make_fdr_facet_labeller <- function(
    annotation
) {
  
  label_strings <- stats::setNames(
    paste0(
      "atop(bold(",
      annotation$gene,
      "), atop(phantom(x), plain(\"",
      annotation$label,
      "\")))"
    ),
    annotation$gene
  )
  
  ggplot2::as_labeller(
    label_strings,
    default = ggplot2::label_parsed
  )
}




###############################################################################
# Panel A: tumor versus adjacent normal expression on log2(TPM + 1)
###############################################################################

tcga_validation_samples <- colnames(
  tpm_all
)[
  substr(
    colnames(
      tpm_all
    ),
    14,
    15
  ) %in%
    c(
      "01",
      "11"
    )
]

tcga_tpm_validation_long <- dplyr::bind_rows(
  lapply(
    names(
      tpm_rows
    ),
    function(gene_name) {
      
      row_id <- unname(
        tpm_rows[
          gene_name
        ]
      )
      
      if (
        length(row_id) != 1L ||
        is.na(row_id) ||
        !row_id %in%
        rownames(
          tpm_all
        )
      ) {
        stop(
          "Could not identify a unique TPM row for ",
          gene_name,
          "."
        )
      }
      
      tibble::tibble(
        sample =
          tcga_validation_samples,
        case_id =
          substr(
            tcga_validation_samples,
            1,
            12
          ),
        gene =
          gene_name,
        condition_code =
          substr(
            tcga_validation_samples,
            14,
            15
          ),
        expression =
          log2(
            as.numeric(
              tpm_all[
                row_id,
                tcga_validation_samples
              ]
            ) + 1
          )
      )
    }
  )
) |>
  dplyr::mutate(
    condition = factor(
      dplyr::if_else(
        condition_code == "11",
        "Normal",
        "Tumor"
      ),
      levels = c(
        "Normal",
        "Tumor"
      )
    )
  ) |>
  dplyr::group_by(
    case_id,
    gene,
    condition
  ) |>
  dplyr::summarise(
    expression = mean(
      expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(
    is.finite(
      expression
    )
  )


# Exact sample sizes represented in Figure 7A after patient-level
# case/condition aggregation and finite-expression filtering.
tcga_case_project_map <- tcga_expression_sample_audit |>
  dplyr::distinct(
    case_id,
    project
  )

tcga_figure7a_case_counts <- tcga_tpm_validation_long |>
  dplyr::left_join(
    tcga_case_project_map,
    by = "case_id"
  ) |>
  dplyr::count(
    gene,
    project,
    condition,
    name = "n_cases"
  ) |>
  dplyr::arrange(
    gene,
    project,
    condition
  )

write_tsv_safe(
  tcga_figure7a_case_counts,
  "TCGA_Figure7A_case_counts.tsv"
)

tcga_figure7a_total_counts <- tcga_tpm_validation_long |>
  dplyr::count(
    gene,
    condition,
    name = "n_cases"
  ) |>
  dplyr::arrange(
    gene,
    condition
  )

write_tsv_safe(
  tcga_figure7a_total_counts,
  "TCGA_Figure7A_total_case_counts.tsv"
)

tcga_figure7a_pair_audit <- tcga_tpm_validation_long |>
  dplyr::distinct(
    case_id,
    gene,
    condition
  ) |>
  dplyr::mutate(
    present = TRUE
  ) |>
  tidyr::pivot_wider(
    names_from = condition,
    values_from = present,
    values_fill = FALSE
  ) |>
  dplyr::group_by(gene) |>
  dplyr::summarise(
    n_cases_with_tumor = sum(Tumor),
    n_cases_with_normal = sum(Normal),
    n_cases_with_both = sum(Tumor & Normal),
    .groups = "drop"
  )

write_tsv_safe(
  tcga_figure7a_pair_audit,
  "TCGA_Figure7A_pair_audit.tsv"
)

message(
  "Figure 7A exact case counts written to TCGA_Figure7A_total_case_counts.tsv ",
  "and TCGA_Figure7A_case_counts.tsv."
)


tcga_validation_stats <- comparison_annotation(
  tcga_tpm_validation_long,
  "condition"
)

if (
  exists(
    "write_tsv_safe",
    inherits = TRUE
  )
) {
  write_tsv_safe(
    tcga_validation_stats,
    "TCGA_focal_gene_expression_tumor_vs_normal_TPM_comparison.tsv"
  )
}


p_tcga_expr_pub <- ggplot2::ggplot(
  tcga_tpm_validation_long,
  ggplot2::aes(
    x = condition,
    y = expression,
    fill = condition
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.58,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.88,
    linewidth = 0.4
  ) +
  ggplot2::facet_wrap(
    ~gene,
    scales = "free_y",
    labeller = make_fdr_facet_labeller(
      tcga_validation_stats
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `Normal` =
        unname(
          PLOT_COLORS[
            "normal"
          ]
        ),
      `Tumor` =
        unname(
          PLOT_COLORS[
            "tumor"
          ]
        )
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.03,
        0.05
      )
    )
  ) +
  ggplot2::labs(
    title =
      "TCGA-COAD/READ expression validation",
    x =
      NULL,
    y =
      expression(
        "log"[2] * "(TPM + 1)"
      ),
    fill =
      "Tissue"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    strip.background =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        face = "plain",
        size = 9.3,
        lineheight = 1.02,
        margin = ggplot2::margin(
          t = 1,
          b = 5
        )
      ),
    axis.text.x =
      ggplot2::element_text(
        size = 8.5
      ),
    legend.position =
      "bottom",
    legend.title =
      ggplot2::element_text(
        face = "bold"
      )
  )



###############################################################################
# Panel B: tumor expression by stage group
###############################################################################

tcga_stage_expression_data_pub <- dplyr::inner_join(
  clinical_stage_data,
  tpm_expression_patient,
  by = "case_id"
) |>
  dplyr::mutate(
    stage_group = factor(
      dplyr::if_else(
        as.character(
          stage_major
        ) %in%
          c(
            "I",
            "II"
          ),
        "Stage I-II",
        "Stage III-IV"
      ),
      levels = c(
        "Stage I-II",
        "Stage III-IV"
      )
    )
  )


stage_expression_long_pub <- tcga_stage_expression_data_pub |>
  dplyr::select(
    case_id,
    stage_group,
    dplyr::all_of(
      FOCAL_GENES
    )
  ) |>
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        FOCAL_GENES
      ),
    names_to =
      "gene",
    values_to =
      "expression"
  ) |>
  dplyr::filter(
    is.finite(
      expression
    )
  )


stage_expression_stats_pub <- comparison_annotation(
  stage_expression_long_pub,
  "stage_group"
)

if (
  exists(
    "write_tsv_safe",
    inherits = TRUE
  )
) {
  write_tsv_safe(
    stage_expression_stats_pub,
    "TCGA_focal_gene_expression_stage_group_comparison.tsv"
  )
}


p_stage_expr_pub <- ggplot2::ggplot(
  stage_expression_long_pub,
  ggplot2::aes(
    x = stage_group,
    y = expression,
    fill = stage_group
  )
) +
  ggplot2::geom_violin(
    trim = FALSE,
    alpha = 0.58,
    linewidth = 0.35
  ) +
  ggplot2::geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.88,
    linewidth = 0.4
  ) +
  ggplot2::facet_wrap(
    ~gene,
    scales = "free_y",
    labeller = make_fdr_facet_labeller(
      stage_expression_stats_pub
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      `Stage I-II` =
        unname(
          PLOT_COLORS[
            "early"
          ]
        ),
      `Stage III-IV` =
        unname(
          PLOT_COLORS[
            "advanced"
          ]
        )
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(
        0.03,
        0.05
      )
    )
  ) +
  ggplot2::labs(
    title =
      "Tumor expression by pathological stage group",
    x =
      NULL,
    y =
      expression(
        "log"[2] * "(TPM + 1)"
      ),
    fill =
      "Stage group"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    strip.background =
      ggplot2::element_blank(),
    strip.text =
      ggplot2::element_text(
        face = "plain",
        size = 9.3,
        lineheight = 1.02,
        margin = ggplot2::margin(
          t = 1,
          b = 5
        )
      ),
    axis.text.x =
      ggplot2::element_text(
        size = 8.5
      ),
    legend.position =
      "bottom",
    legend.title =
      ggplot2::element_text(
        face = "bold"
      )
  )



###############################################################################
# Panel C: stage-dependent Cox interaction forest plot
###############################################################################

stage_forest_df_pub <- dplyr::bind_rows(
  stage_interaction_summary |>
    dplyr::transmute(
      gene =
        gene,
      stage_group =
        "Stage I-II",
      HR =
        HR_Stage_I_II,
      CI_low =
        CI_low_Stage_I_II,
      CI_high =
        CI_high_Stage_I_II,
      interaction_FDR =
        interaction_FDR
    ),
  stage_interaction_summary |>
    dplyr::transmute(
      gene =
        gene,
      stage_group =
        "Stage III-IV",
      HR =
        HR_Stage_III_IV,
      CI_low =
        CI_low_Stage_III_IV,
      CI_high =
        CI_high_Stage_III_IV,
      interaction_FDR =
        interaction_FDR
    )
) |>
  dplyr::mutate(
    gene =
      factor(
        gene,
        levels =
          rev(
            FOCAL_GENES
          )
      ),
    stage_group =
      factor(
        stage_group,
        levels = c(
          "Stage I-II",
          "Stage III-IV"
        )
      )
  )


forest_x_min <- min(
  stage_forest_df_pub$CI_low,
  na.rm = TRUE
)

forest_ci_max <- max(
  stage_forest_df_pub$CI_high,
  na.rm = TRUE
)

forest_x_max <- forest_ci_max * 1.48
forest_label_x <- forest_x_max / 1.035


stage_fdr_labels_pub <- stage_interaction_summary |>
  dplyr::transmute(
    gene =
      factor(
        gene,
        levels =
          rev(
            FOCAL_GENES
          )
      ),
    label =
      paste0(
        "Interaction FDR = ",
        vapply(
          interaction_FDR,
          format_fdr,
          character(1)
        )
      ),
    label_x =
      forest_label_x
  )


p_stage_forest_pub <- ggplot2::ggplot(
  stage_forest_df_pub,
  ggplot2::aes(
    x = HR,
    y = gene,
    color = stage_group
  )
) +
  ggplot2::geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    height = 0.13,
    linewidth = 0.8,
    position =
      ggplot2::position_dodge(
        width = 0.42
      )
  ) +
  ggplot2::geom_point(
    size = 3,
    position =
      ggplot2::position_dodge(
        width = 0.42
      )
  ) +
  ggplot2::geom_text(
    data =
      stage_fdr_labels_pub,
    ggplot2::aes(
      x = label_x,
      y = gene,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    fontface = "bold",
    size = 3.05,
    color = "grey20"
  ) +
  ggplot2::scale_color_manual(
    values = c(
      `Stage I-II` =
        unname(
          PLOT_COLORS[
            "early"
          ]
        ),
      `Stage III-IV` =
        unname(
          PLOT_COLORS[
            "advanced"
          ]
        )
    )
  ) +
  ggplot2::scale_x_log10(
    limits = c(
      forest_x_min * 0.90,
      forest_x_max
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.01,
        0
      )
    )
  ) +
  ggplot2::labs(
    title =
      "Stage-dependent overall-survival associations",
    x =
      "Hazard ratio per 1-SD increase in log2(TPM+1)",
    y =
      NULL,
    color =
      "Stage group"
  ) +
  ggplot2::theme_classic(
    base_size = 10.5
  ) +
  ggplot2::theme(
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 12
      ),
    axis.title.x =
      ggplot2::element_text(
        face = "bold"
      ),
    axis.text.y =
      ggplot2::element_text(
        face = "bold"
      ),
    legend.position =
      "bottom",
    legend.title =
      ggplot2::element_text(
        face = "bold"
      ),
    plot.margin =
      ggplot2::margin(
        6,
        8,
        5,
        6
      )
  )


###############################################################################
# Figure assembly
###############################################################################

top_row <- p_tcga_expr_pub + p_stage_expr_pub +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )


# The forest plot is centered and narrower than the full top row.
bottom_row <- patchwork::plot_spacer() +
  p_stage_forest_pub +
  patchwork::plot_spacer() +
  patchwork::plot_layout(
    widths = c(
      0.17,
      1.66,
      0.17
    )
  )


figure7_tcga_pub <- top_row / bottom_row +
  patchwork::plot_layout(
    heights = c(
      1,
      0.82
    )
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.tag =
        ggplot2::element_text(
          face = "bold",
          size = 13
        ),
      plot.tag.position =
        "topleft"
    )
  )



# Figure 8 palette
km_low_colour <- "#0072B2"
km_high_colour <- "#D7191C"
km_source <- km_plots
panel_info <- data.frame(
  key = c(
    "Stage_I_II__AQP8",
    "Stage_III_IV__AQP8",
    "Stage_I_II__GUCA2A",
    "Stage_III_IV__GUCA2A",
    "Stage_I_II__MS4A12",
    "Stage_III_IV__MS4A12"
  ),
  gene = c(
    "AQP8",
    "AQP8",
    "GUCA2A",
    "GUCA2A",
    "MS4A12",
    "MS4A12"
  ),
  stage = c(
    "Stage I-II",
    "Stage III-IV",
    "Stage I-II",
    "Stage III-IV",
    "Stage I-II",
    "Stage III-IV"
  ),
  tag = c(
    "A",
    "D",
    "B",
    "E",
    "C",
    "F"
  ),
  stringsAsFactors = FALSE
)

panel_info$p_value <- vapply(
  seq_len(
    nrow(panel_info)
  ),
  function(i) {
    
    hit <- km_summary[
      km_summary$gene == panel_info$gene[i] &
        km_summary$stage_group == panel_info$stage[i],
      ,
      drop = FALSE
    ]
    
    if (nrow(hit) != 1L) {
      return(NA_real_)
    }
    
    as.numeric(
      hit$logrank_p[1]
    )
  },
  numeric(1)
)


missing_keys <- setdiff(
  panel_info$key,
  names(km_source)
)

if (length(missing_keys) > 0) {
  stop(
    "Missing KM panel(s): ",
    paste(
      missing_keys,
      collapse = ", "
    )
  )
}


style_km_panel <- function(
    km_object,
    gene,
    stage,
    tag,
    p_value
) {
  
  km_object$plot <-
    km_object$plot +
    ggplot2::scale_color_manual(
      values = c(
        km_low_colour,
        km_high_colour
      )
    ) +
    ggplot2::labs(
      title = gene,
      subtitle = stage,
      tag = tag
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 2
        )
      ),
      
      plot.subtitle = ggplot2::element_text(
        face = "plain",
        size = 9.6,
        color = "grey30",
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 7
        )
      ),
      
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      
      plot.tag.position = c(
        0.01,
        1.02
      ),
      
      plot.margin = ggplot2::margin(
        t = 12,
        r = 22,
        b = 7,
        l = 22
      ),
      
      legend.position = "top",
      legend.margin = ggplot2::margin(
        b = 4
      )
    )
  
  
  if (!is.null(km_object$table)) {
    
    km_object$table <-
      km_object$table +
      ggplot2::scale_color_manual(
        values = c(
          km_low_colour,
          km_high_colour
        )
      ) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(
          t = 3,
          r = 22,
          b = 13,
          l = 22
        ),
        
        plot.title = ggplot2::element_text(
          face = "plain",
          size = 10,
          margin = ggplot2::margin(
            b = 4
          )
        )
      )
  }
  
  
  # Bold only statistically significant log-rank p-values.
  # survminer adds the p-value as a text layer to the KM plot.
  if (
    is.finite(p_value) &&
    length(km_object$plot$layers) > 0
  ) {
    
    p_fontface <- if (
      p_value < 0.05
    ) {
      "bold"
    } else {
      "plain"
    }
    
    for (
      layer_index in
      seq_along(
        km_object$plot$layers
      )
    ) {
      
      layer <- km_object$plot$layers[[layer_index]]
      
      if (
        inherits(
          layer$geom,
          "GeomText"
        )
      ) {
        km_object$plot$layers[[layer_index]]$aes_params$fontface <-
          p_fontface
      }
    }
  }
  
  
  km_object
}


km_publication <- vector(
  "list",
  nrow(
    panel_info
  )
)

for (i in seq_len(nrow(panel_info))) {
  
  km_publication[[i]] <-
    style_km_panel(
      km_object =
        km_source[[panel_info$key[i]]],
      gene =
        panel_info$gene[i],
      stage =
        panel_info$stage[i],
      tag =
        panel_info$tag[i],
      p_value =
        panel_info$p_value[i]
    )
}


# arrange_ggsurvplots fills the grid column-wise.
# This order therefore produces:
# A B C
# D E F
figure8_pub <-
  survminer::arrange_ggsurvplots(
    km_publication,
    print = FALSE,
    ncol = 3,
    nrow = 2
  )




###############################################################################
# Save publication figures
###############################################################################

save_plot(
  figure1_pub,
  "Figure_1_Meta_DEG",
  width = 11,
  height = 5,
  dpi = 400
)

save_manuscript_plot(
  figure1_pub,
  "Figure_1_Meta_DEG",
  width = 11,
  height = 5,
  dpi = 600
)

save_plot(
  figure2_pub,
  "Figure_2_GO_enrichment",
  width = 12.5,
  height = 12.2,
  dpi = 400
)

save_manuscript_plot(
  figure2_pub,
  "Figure_2_GO_enrichment",
  width = 12.5,
  height = 12.2,
  dpi = 600
)

save_plot(
  figure3_pub,
  "Figure_3_KEGG_enrichment",
  width = 12,
  height = 6.8,
  dpi = 400
)

save_manuscript_plot(
  figure3_pub,
  "Figure_3_KEGG_enrichment",
  width = 12,
  height = 6.8,
  dpi = 600
)

save_plot(
  figure4_pub,
  "Figure_4_PPI_MCODE_modules",
  width = 11.5,
  height = 13.2,
  dpi = 400
)

save_manuscript_plot(
  figure4_pub,
  "Figure_4_PPI_MCODE_modules",
  width = 11.5,
  height = 13.2,
  dpi = 600
)

save_plot(
  p_sft1_pub +
    p_sft2_pub +
    patchwork::plot_layout(
      widths = c(1, 1)
    ) +
    patchwork::plot_annotation(
      tag_levels = "A"
    ),
  "WGCNA_soft_threshold_diagnostics",
  width = 10,
  height = 5,
  dpi = 400
)

save_plot(
  p_modtrait_pub,
  "WGCNA_module_trait_heatmap_paired",
  width = 6.2,
  height = 6.0,
  dpi = 400
)

save_plot(
  p_modtrait_pub,
  "Figure_5D_WGCNA_module_trait",
  width = 6.2,
  height = 6.0,
  dpi = 400
)

save_plot(
  figure5_wgcna_pub,
  "Figure_5_WGCNA",
  width = 11.5,
  height = 13.2,
  dpi = 400
)

save_manuscript_plot(
  figure5_wgcna_pub,
  "Figure_5_WGCNA",
  width = 11.5,
  height = 13.2,
  dpi = 600
)

save_plot(
  figure6_wgcna_pub,
  "Figure_6_WGCNA_hubs_and_overlap",
  width = 12,
  height = 9.4,
  dpi = 400
)

save_manuscript_plot(
  figure6_wgcna_pub,
  "Figure_6_WGCNA_hubs_and_overlap",
  width = 12,
  height = 9.4,
  dpi = 600
)

save_plot(
  p_tcga_expr_pub,
  "TCGA_focal_gene_expression",
  width = 8.5,
  height = 4.7,
  dpi = 400
)

save_plot(
  p_stage_expr_pub,
  "TCGA_focal_gene_expression_by_stage_group",
  width = 8.5,
  height = 4.7,
  dpi = 400
)

save_plot(
  p_stage_forest_pub,
  "TCGA_TPM_robust_stage_forest",
  width = 7.8,
  height = 4.7,
  dpi = 400
)

save_plot(
  figure7_tcga_pub,
  "Figure_7_TCGA_validation_and_survival",
  width = 12,
  height = 9.1,
  dpi = 400
)

save_manuscript_plot(
  figure7_tcga_pub,
  "Figure_7_TCGA_validation_and_survival",
  width = 12,
  height = 9.1,
  dpi = 600
)

save_survival_publication <- function(
    plot,
    directory,
    filename,
    width = 16.5,
    height = 11.2,
    dpi = 400
) {
  
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  grDevices::cairo_pdf(
    filename = file.path(
      directory,
      paste0(
        filename,
        ".pdf"
      )
    ),
    width = width,
    height = height
  )
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  grDevices::png(
    filename = file.path(
      directory,
      paste0(
        filename,
        ".png"
      )
    ),
    width = width * dpi,
    height = height * dpi,
    res = dpi
  )
  
  print(
    plot
  )
  
  grDevices::dev.off()
}

save_survival_publication(
  figure8_pub,
  DIR_FIG,
  "Figure_8_TCGA_stage_survival",
  dpi = 400
)

save_survival_publication(
  figure8_pub,
  DIR_MANUSCRIPT,
  "Figure_8_TCGA_stage_survival",
  dpi = 600
)




###############################################################################
# 29. SUPPLEMENTARY FIGURE HARMONIZATION
###############################################################################

# Supplementary enrichment figures use the same visual grammar as Figures 2–3:
# point size = gene count, color intensity = -log10(FDR), blue for downregulated
# or normal-associated sets, and red/turquoise for the corresponding comparison.

make_supp_enrichment_panel <- function(
    data,
    gene_set_name,
    tag,
    title_text,
    low_colour,
    high_colour,
    size_limits
) {
  
  d <- data |>
    dplyr::filter(
      gene_set == gene_set_name,
      is.finite(gene_ratio),
      is.finite(count),
      is.finite(FDR),
      FDR > 0
    ) |>
    dplyr::mutate(
      neglog10_fdr = -log10(FDR)
    )
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(
      x = gene_ratio,
      y = stats::reorder(
        term,
        gene_ratio
      ),
      size = count,
      fill = neglog10_fdr
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      color = "grey25",
      stroke = 0.24,
      alpha = 0.97
    ) +
    ggplot2::facet_grid(
      collection ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_fill_gradient(
      low = low_colour,
      high = high_colour,
      name = expression(
        -log[10](FDR)
      )
    ) +
    ggplot2::scale_size_continuous(
      name = "Genes",
      range = c(
        2.4,
        7.0
      ),
      limits = size_limits,
      breaks = scales::pretty_breaks(
        n = 3
      )
    ) +
    ggplot2::labs(
      title = title_text,
      tag = tag,
      x = "Gene ratio",
      y = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 10.5
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 11.5,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        face = "bold",
        size = 13
      ),
      plot.tag.position = c(
        0.01,
        1.015
      ),
      strip.background = ggplot2::element_blank(),
      strip.text.y.right = ggplot2::element_text(
        face = "bold",
        angle = 0,
        size = 8.5,
        margin = ggplot2::margin(
          l = 5
        )
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 6.8
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 8.5
      ),
      legend.text = ggplot2::element_text(
        size = 7.5
      ),
      plot.margin = ggplot2::margin(
        8,
        8,
        6,
        8
      )
    )
}


if (
  exists(
    "deg_enrichment_supp",
    inherits = TRUE
  ) &&
  nrow(
    deg_enrichment_supp
  ) > 0
) {
  
  deg_supp_size_limits <- range(
    deg_enrichment_supp$count,
    na.rm = TRUE
  )
  
  p_deg_down_supp <- make_supp_enrichment_panel(
    deg_enrichment_supp,
    "Downregulated DEGs",
    "A",
    "Downregulated DEGs",
    "#DCEAF7",
    "#2166AC",
    deg_supp_size_limits
  )
  
  p_deg_up_supp <- make_supp_enrichment_panel(
    deg_enrichment_supp,
    "Upregulated DEGs",
    "B",
    "Upregulated DEGs",
    "#F7D9D5",
    "#B2182B",
    deg_supp_size_limits
  )
  
  figure_s07_pub <-
    p_deg_down_supp +
    p_deg_up_supp +
    patchwork::plot_layout(
      widths = c(
        1,
        1
      )
    ) +
    patchwork::plot_annotation(
      title = "Functional enrichment of meta-analysis DEGs",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 13,
          hjust = 0.5
        )
      )
    )
  
  save_supplementary_plot(
    figure_s07_pub,
    "Figure_S07_DEG_full_enrichment",
    width = 14,
    height = 18,
    dpi = 600
  )
}


if (
  exists(
    "wgcna_enrichment_supp",
    inherits = TRUE
  ) &&
  nrow(
    wgcna_enrichment_supp
  ) > 0
) {
  
  wgcna_supp_size_limits <- range(
    wgcna_enrichment_supp$count,
    na.rm = TRUE
  )
  
  turquoise_set_name <- paste0(
    turquoise_module,
    " module"
  )
  
  blue_set_name <- paste0(
    blue_module,
    " module"
  )
  
  p_wgcna_turquoise_supp <- make_supp_enrichment_panel(
    wgcna_enrichment_supp,
    turquoise_set_name,
    "A",
    paste0(
      "Module ",
      turquoise_module
    ),
    "#D9F1F0",
    "#008C95",
    wgcna_supp_size_limits
  )
  
  p_wgcna_blue_supp <- make_supp_enrichment_panel(
    wgcna_enrichment_supp,
    blue_set_name,
    "B",
    paste0(
      "Module ",
      blue_module
    ),
    "#DCEAF7",
    "#2166AC",
    wgcna_supp_size_limits
  )
  
  figure_s08_pub <-
    p_wgcna_turquoise_supp +
    p_wgcna_blue_supp +
    patchwork::plot_layout(
      widths = c(
        1,
        1
      )
    ) +
    patchwork::plot_annotation(
      title = "Functional enrichment of WGCNA modules",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          size = 13,
          hjust = 0.5
        )
      )
    )
  
  save_supplementary_plot(
    figure_s08_pub,
    "Figure_S08_WGCNA_module_enrichment",
    width = 14,
    height = 18,
    dpi = 600
  )
}


if (
  exists(
    "ppi_top",
    inherits = TRUE
  ) &&
  nrow(
    ppi_top
  ) > 0
) {
  
  p_ppi_centrality_pub <- ggplot2::ggplot(
    ppi_top,
    ggplot2::aes(
      x = final_score,
      y = stats::reorder(
        symbol,
        final_score
      ),
      color = focal
    )
  ) +
    ggplot2::geom_point(
      size = 2.9,
      alpha = 0.95
    ) +
    ggplot2::scale_color_manual(
      values = c(
        `FALSE` = "grey68",
        `TRUE` = "#B2182B"
      ),
      labels = c(
        `FALSE` = "Other gene",
        `TRUE` = "Focal gene"
      )
    ) +
    ggplot2::labs(
      title = "Highest-ranked STRING/PPI genes",
      x = "Final Centrality score",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_classic(
      base_size = 10.5
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      legend.position = "bottom",
      axis.text.y = ggplot2::element_text(
        size = 8
      )
    )
  
  save_supplementary_plot(
    p_ppi_centrality_pub,
    "Figure_S09_PPI_centrality",
    width = 7.5,
    height = 8,
    dpi = 600
  )
}


if (
  exists(
    "p_overall_cox",
    inherits = TRUE
  )
) {
  
  p_overall_cox_pub <-
    p_overall_cox +
    ggplot2::labs(
      subtitle = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 12
      ),
      axis.title.x = ggplot2::element_text(
        face = "bold"
      )
    )
  
  save_supplementary_plot(
    p_overall_cox_pub,
    "Figure_S11_TCGA_overall_Cox",
    width = 7,
    height = 4.5,
    dpi = 600
  )
}




###############################################################################
# 30. FIGURE MANIFEST AND REPRODUCIBILITY
###############################################################################

analysis_parameters <- tibble::tibble(
  parameter = c(
    "Meta-analysis method",
    "Meta-analysis FDR cutoff",
    "Meta-analysis absolute log2FC cutoff",
    "Minimum contributing GEO studies",
    "WGCNA network type",
    "WGCNA TOM type",
    "WGCNA selected soft-threshold power",
    "WGCNA scale-free R2 target",
    "WGCNA minimum module size",
    "WGCNA merge cut height",
    "WGCNA hub |MM| cutoff",
    "WGCNA hub GS cutoff",
    "STRING version",
    "STRING score threshold",
    "TCGA focal-gene expression scale",
    "Survival primary expression model",
    "Survival stage-interaction comparison"
  ),
  value = c(
    META_METHOD,
    as.character(META_FDR_CUTOFF),
    as.character(META_LOGFC_CUTOFF),
    as.character(META_MIN_STUDIES),
    "signed",
    "signed",
    as.character(WGCNA_SOFT_POWER_USED),
    as.character(WGCNA_SFT_R2_TARGET),
    as.character(WGCNA_MIN_MODULE_SIZE),
    as.character(WGCNA_MERGE_CUT),
    as.character(WGCNA_MM_CUTOFF),
    as.character(WGCNA_GS_CUTOFF),
    STRING_VERSION,
    as.character(STRING_SCORE_THRESHOLD),
    "log2(TPM + 1)",
    "Continuous expression per 1-SD increase; adjusted for age, sex, project, and exact pathological stage",
    "Likelihood-ratio test of expression-by-stage interaction; Stage I-II versus Stage III-IV"
  )
)

write_tsv_safe(
  analysis_parameters,
  "Analysis_parameters.tsv"
)


figure_manifest <- dplyr::bind_rows(
  tibble::tibble(
    folder = "Manuscript",
    file = list.files(
      DIR_MANUSCRIPT,
      full.names = FALSE
    )
  ),
  tibble::tibble(
    folder = "Supplementary",
    file = list.files(
      DIR_SUPPLEMENTARY,
      full.names = FALSE
    )
  )
) |>
  dplyr::arrange(
    folder,
    file
  )

write_tsv_safe(
  figure_manifest,
  "Figure_manifest.tsv"
)

capture.output(
  sessionInfo(),
  file = file.path(
    DIR_RESULTS,
    "sessionInfo.txt"
  )
)

message(
  "\n============================================================"
)

message(
  "Analysis complete"
)

message(
  "Results: ",
  DIR_RESULTS
)

message(
  "Tables: ",
  DIR_TABLE
)

message(
  "Figures: ",
  DIR_FIG
)

message(
  "Manuscript figures: ",
  DIR_MANUSCRIPT
)

message(
  "Supplementary figures: ",
  DIR_SUPPLEMENTARY
)

message(
  "============================================================\n"
)
