#!/usr/bin/env Rscript

# ============================================================================
# COMPASS baseline benchmark
# FINAL THESIS = canonical COMPASS source of truth
#
# COMPASS score:
#   taken directly from data/human_normal_all.csv (same file used by thesis ROC)
#
# Baselines:
#   reconstructed from the final RNA and ChIP modality-level scores,
#   without changing upstream RNA/ChIP processing.
#
# This script first reproduces the original thesis ROC exactly. If that sanity
# check fails materially, baseline comparisons should NOT be interpreted yet.
# ============================================================================

source("config.R")

required_pkgs <- c("data.table", "dplyr", "tidyr", "readr", "pROC", "ggplot2")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nRun:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(pROC)
  library(ggplot2)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

log_lines <- character()
say <- function(...) {
  x <- paste0(...)
  message(x)
  log_lines <<- c(log_lines, x)
}

die_if_missing <- function(path, label = basename(path)) {
  if (!file.exists(path)) stop("Missing required file: ", label, "\nExpected: ", path)
}

normalize_gene <- function(x) toupper(trimws(as.character(x)))

normalize_species <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z[z %in% c("hs", "hg38", "homo sapiens", "homo_sapiens")] <- "human"
  z[z %in% c("mm", "mm10", "mm39", "mus musculus", "mus_musculus")] <- "mouse"
  z
}

normalize_status <- function(x) tolower(trimws(as.character(x)))

normalize_tissue <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("_", "-", z, fixed = TRUE)
  z <- gsub("[[:space:]]+", "-", z)
  z <- gsub("-+", "-", z)
  z[z == "testis"] <- "testicle"
  z[z == "head-and neck"] <- "head-and-neck"
  z
}

load_rda_list <- function(path) {
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  mget(nm, envir = e, inherits = FALSE)
}

require_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(label, " missing columns: ", paste(miss, collapse = ", "))
  }
}

read_gmt_genes <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nchar(ln) > 0][1]
  normalize_gene(strsplit(ln, "\t")[[1]][-(1:2)])
}

rna_conf_score <- function(x) {
  dplyr::case_when(
    tolower(as.character(x)) == "high" ~ 1.0,
    tolower(as.character(x)) == "medium" ~ 0.4,
    tolower(as.character(x)) == "low" ~ 0.2,
    TRUE ~ 0
  )
}

chip_bind_score <- function(x) {
  dplyr::case_when(
    as.character(x) == "High confidence binding" ~ 1.0,
    as.character(x) == "Medium-high confidence" ~ 0.8,
    as.character(x) == "Medium confidence" ~ 0.4,
    as.character(x) == "Low-medium confidence" ~ 0.3,
    as.character(x) == "Low confidence" ~ 0.15,
    TRUE ~ 0
  )
}

auc_ci_thesis <- function(labels, scores, n = THESIS_ROC_BOOT_N) {
  r <- pROC::roc(
    labels, scores,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  ci <- as.numeric(
    pROC::ci.auc(
      r,
      method = "bootstrap",
      boot.n = n,
      boot.stratified = FALSE,
      progress = "none"
    )
  )
  list(
    roc = r,
    auc = as.numeric(pROC::auc(r)),
    lo = ci[1],
    hi = ci[3]
  )
}

average_precision_grouped <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- score[ok]
  npos <- sum(y == 1)
  if (npos == 0) return(NA_real_)

  o <- order(score, decreasing = TRUE)
  y <- y[o]
  score <- score[o]

  grp <- cumsum(c(TRUE, diff(score) != 0))
  pos_by <- as.numeric(tapply(y, grp, sum))
  n_by <- as.numeric(tapply(y, grp, length))

  cum_pos <- cumsum(pos_by)
  cum_n <- cumsum(n_by)
  precision <- cum_pos / cum_n
  recall <- cum_pos / npos
  delta_recall <- c(recall[1], diff(recall))

  sum(delta_recall * precision)
}

topk_metrics <- function(gene, y, score, k) {
  kk <- min(as.integer(k), length(score))
  # deterministic tie-breaking only for Top-K extraction
  o <- order(-score, gene)
  idx <- o[seq_len(kk)]
  tp <- sum(y[idx] == 1)
  prevalence <- mean(y == 1)
  precision <- tp / kk
  recall <- tp / sum(y == 1)
  enrichment <- precision / prevalence
  tibble(
    k = kk,
    true_positives = tp,
    precision = precision,
    recall = recall,
    enrichment = enrichment
  )
}

# ============================================================================
# 0. Input checks
# ============================================================================
die_if_missing(COMPASS_CSV)
for (x in GMT) die_if_missing(x)
die_if_missing(RNA_RDA)
die_if_missing(CHIP_RDA)

say("RUN_MODE = ", RUN_MODE)
say("N_BOOT = ", N_BOOT)
say("Canonical COMPASS file = ", COMPASS_CSV)

# ============================================================================
# 1. Reproduce the original thesis ROC exactly
# ============================================================================
say("\n=== STEP 1: reproduce thesis ROC ===")

dt <- data.table::fread(COMPASS_CSV)
require_cols(dt, c("Gene", "Gene_biotype", "Tissue", "PTS"), "human_normal_all.csv")

dt <- dt[Gene_biotype == GENE_BIOTYPE]
dt[, Gene := normalize_gene(Gene)]

g_all <- dt[, .(PTS = max(PTS)), by = Gene]
g_ne  <- dt[Tissue != EXCLUDE_TISSUE, .(PTS = max(PTS)), by = Gene]
sets <- lapply(GMT, read_gmt_genes)

roc_rows <- lapply(names(sets), function(k) {
  s <- sets[[k]]
  lab <- as.integer(g_ne$Gene %in% s)
  rc <- auc_ci_thesis(lab, g_ne$PTS)

  lab_all <- as.integer(g_all$Gene %in% s)
  auc_all <- as.numeric(
    pROC::auc(
      pROC::roc(
        lab_all, g_all$PTS,
        levels = c(0,1),
        direction = "<",
        quiet = TRUE
      )
    )
  )

  data.table(
    set = k,
    npos = sum(lab),
    n = length(lab),
    auc = rc$auc,
    lo = rc$lo,
    hi = rc$hi,
    auc_all = auc_all
  )
})

roc_repro <- data.table::rbindlist(roc_rows)
data.table::fwrite(
  roc_repro,
  file.path(OUT_DIR, "01_thesis_ROC_reproduction.tsv"),
  sep = "\t"
)

say("Protein-coding gene universe (non-ESC aggregate) = ", nrow(g_ne))
for (i in seq_len(nrow(roc_repro))) {
  say(
    roc_repro$set[i],
    ": positives=", roc_repro$npos[i],
    ", AUC=", round(roc_repro$auc[i], 4),
    ", 95% CI=", round(roc_repro$lo[i], 4), "–", round(roc_repro$hi[i], 4),
    ", AUC_all=", round(roc_repro$auc_all[i], 4)
  )
}

sanity <- tibble(
  check = c(
    "universe_n",
    paste0("positives_", names(EXPECTED_POSITIVES)),
    paste0("AUC_", names(EXPECTED_AUC))
  ),
  expected = c(
    EXPECTED_N,
    as.numeric(EXPECTED_POSITIVES),
    as.numeric(EXPECTED_AUC)
  ),
  observed = c(
    nrow(g_ne),
    vapply(names(EXPECTED_POSITIVES), function(k) {
      roc_repro$npos[roc_repro$set == k][1]
    }, numeric(1)),
    vapply(names(EXPECTED_AUC), function(k) {
      roc_repro$auc[roc_repro$set == k][1]
    }, numeric(1))
  )
) %>%
  mutate(
    difference = observed - expected,
    pass = case_when(
      check == "universe_n" ~ observed == expected,
      grepl("^positives_", check) ~ observed == expected,
      grepl("^AUC_", check) ~ abs(difference) <= AUC_TOLERANCE,
      TRUE ~ FALSE
    )
  )

readr::write_tsv(sanity, file.path(OUT_DIR, "02_sanity_checks.tsv"))

if (!all(sanity$pass)) {
  say(
    "WARNING: thesis ROC sanity check is not fully reproduced. ",
    "Baseline outputs will still be generated for debugging, but should not be ",
    "interpreted for the manuscript until the mismatch is resolved."
  )
} else {
  say("PASS: thesis ROC sanity checks reproduced.")
}


# ============================================================================
# 1b. Final-journal ROC confidence intervals (2000 bootstrap)
#     Canonical thesis reproduction above remains unchanged at 1000 bootstrap.
# ============================================================================
say("\n=== STEP 1b: final-journal ROC confidence intervals ===")

auc_ci_journal <- function(labels, scores, n = JOURNAL_ROC_BOOT_N) {
  r <- pROC::roc(
    labels, scores,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  ci <- as.numeric(
    pROC::ci.auc(
      r,
      method = "bootstrap",
      boot.n = n,
      boot.stratified = FALSE,
      progress = "none"
    )
  )
  list(
    auc = as.numeric(pROC::auc(r)),
    lo = ci[1],
    hi = ci[3]
  )
}

journal_roc_rows <- lapply(names(sets), function(k) {
  s <- sets[[k]]
  lab <- as.integer(g_ne$Gene %in% s)
  rc <- auc_ci_journal(lab, g_ne$PTS, JOURNAL_ROC_BOOT_N)
  data.table(
    set = k,
    npos = sum(lab),
    n = length(lab),
    auc = rc$auc,
    lo = rc$lo,
    hi = rc$hi,
    boot_n = JOURNAL_ROC_BOOT_N
  )
})

journal_roc <- data.table::rbindlist(journal_roc_rows)
data.table::fwrite(
  journal_roc,
  file.path(OUT_DIR, "02b_journal_ROC_CI_2000.tsv"),
  sep = "\t"
)

for (i in seq_len(nrow(journal_roc))) {
  say(
    "Journal CI ", journal_roc$set[i],
    ": AUC=", round(journal_roc$auc[i], 4),
    ", 95% CI=", round(journal_roc$lo[i], 4), "–", round(journal_roc$hi[i], 4),
    ", B=", JOURNAL_ROC_BOOT_N
  )
}


# Canonical universe and canonical COMPASS score
canonical <- as_tibble(g_ne) %>%
  rename(gene = Gene, COMPASS = PTS) %>%
  mutate(gene = normalize_gene(gene)) %>%
  distinct(gene, .keep_all = TRUE)

# ============================================================================
# 2. Load final RNA and ChIP modality-level scores
# ============================================================================
say("\n=== STEP 2: load modality-level final scores ===")

rna_objs <- load_rda_list(RNA_RDA)
chip_objs <- load_rda_list(CHIP_RDA)

rna_expected <- c(
  "human_normal_rna_final_score",
  "human_disease_rna_final_score",
  "mouse_normal_rna_final_score",
  "mouse_disease_rna_final_score"
)

chip_expected <- c(
  "human_normal_chip_final_score",
  "human_disease_chip_final_score",
  "mouse_normal_chip_final_score",
  "mouse_disease_chip_final_score"
)

if (length(setdiff(rna_expected, names(rna_objs))) > 0) {
  stop(
    "RNA RDA is missing expected objects: ",
    paste(setdiff(rna_expected, names(rna_objs)), collapse = ", ")
  )
}
if (length(setdiff(chip_expected, names(chip_objs))) > 0) {
  stop(
    "ChIP RDA is missing expected objects: ",
    paste(setdiff(chip_expected, names(chip_objs)), collapse = ", ")
  )
}

rna0 <- rna_objs[["human_normal_rna_final_score"]]
chip0 <- chip_objs[["human_normal_chip_final_score"]]

require_cols(
  rna0,
  c("tissue", "gene", "final_scaled_score", "confidence_level"),
  "human_normal_rna_final_score"
)
require_cols(
  chip0,
  c("tissue", "gene", "final_scaled_score", "binding_class"),
  "human_normal_chip_final_score"
)

rna <- as_tibble(rna0) %>%
  transmute(
    tissue = normalize_tissue(tissue),
    gene = normalize_gene(gene),
    rna_scaled = as.numeric(final_scaled_score),
    rna_confidence = tolower(as.character(confidence_level)),
    rna_confidence_score = rna_conf_score(rna_confidence),
    rna_component = rna_scaled * rna_confidence_score,
    legacy_rna_regulation = if ("regulation" %in% names(rna0)) {
      tolower(as.character(rna0$regulation))
    } else {
      NA_character_
    }
  )

chip <- as_tibble(chip0) %>%
  transmute(
    tissue = normalize_tissue(tissue),
    gene = normalize_gene(gene),
    chip_scaled = as.numeric(final_scaled_score),
    binding_class = as.character(binding_class),
    chip_binding_score = chip_bind_score(binding_class),
    chip_component = chip_scaled * chip_binding_score
  )

# Exclude ESC for the benchmark, exactly as the original ROC did
rna <- rna %>% filter(tissue != normalize_tissue(EXCLUDE_TISSUE))
chip <- chip %>% filter(tissue != normalize_tissue(EXCLUDE_TISSUE))

# Optional direction recovery
# NOTE:
# The saved 1.3-rna_freq.rda produced by find_common_genes() contains the
# tissue-level direction in column `regulation`, because the function renames
# `dominant_direction` -> `regulation` immediately before returning.
# This reader accepts BOTH names so it is compatible with historical objects.
direction_source <- "legacy_rna_final_score"
if (file.exists(RNA_FREQ_RDA)) {
  say("Optional RNA frequency file found: ", RNA_FREQ_RDA)
  freq_objs <- load_rda_list(RNA_FREQ_RDA)
  target_name <- "human_normal_freq_list"

  if (target_name %in% names(freq_objs)) {
    lst <- freq_objs[[target_name]]

    freq_rows <- lapply(names(lst), function(nm) {
      x <- as_tibble(lst[[nm]])
      names(x) <- tolower(names(x))

      if (!"gene" %in% names(x)) {
        warning("Skipping ", nm, ": no gene column.")
        return(NULL)
      }

      # Current saved object: regulation.
      # Older/intermediate object, if any: dominant_direction.
      direction_col <- if ("regulation" %in% names(x)) {
        "regulation"
      } else if ("dominant_direction" %in% names(x)) {
        "dominant_direction"
      } else {
        warning(
          "Skipping ", nm,
          ": neither regulation nor dominant_direction was found. Columns: ",
          paste(names(x), collapse = ", ")
        )
        return(NULL)
      }

      tissue_from_name <- sub("^human_normal_", "", nm)

      tibble(
        tissue = normalize_tissue(tissue_from_name),
        gene = normalize_gene(x$gene),
        source_direction = tolower(as.character(x[[direction_col]]))
      )
    })

    # Remove NULL entries before bind_rows. This avoids the previous
    # zero-column data frame / distinct(tissue, gene) failure.
    freq_rows <- Filter(Negate(is.null), freq_rows)

    if (length(freq_rows) > 0) {
      freq_dir <- bind_rows(freq_rows) %>%
        filter(!is.na(gene), gene != "") %>%
        distinct(tissue, gene, .keep_all = TRUE) %>%
        mutate(
          rna_direction = case_when(
            source_direction == "up" ~ "up",
            source_direction == "down" ~ "down",
            source_direction %in% c("mixed", "unchanged") ~ "no_dominant_direction",
            TRUE ~ "no_dominant_direction"
          )
        )

      rna <- rna %>%
        left_join(
          freq_dir %>% select(tissue, gene, rna_direction),
          by = c("tissue", "gene")
        )

      direction_source <- "1.3-rna_freq.rda"
      say(
        "Recovered RNA directions from 1.3-rna_freq.rda: ",
        nrow(freq_dir), " gene-tissue rows."
      )
      say(
        "Recovered direction counts: ",
        paste(
          names(table(freq_dir$rna_direction)),
          as.integer(table(freq_dir$rna_direction)),
          collapse = "; "
        )
      )
    } else {
      warning(
        "1.3-rna_freq.rda was found, but no usable gene/direction rows were ",
        "parsed. Falling back to regulation in 1.5-rna_final_score.rda."
      )
    }
  } else {
    warning(
      "1.3-rna_freq.rda does not contain object human_normal_freq_list. ",
      "Falling back to regulation in 1.5-rna_final_score.rda."
    )
  }
}

# Fill any gene-tissue rows that were not found in the optional frequency file.
if (!"rna_direction" %in% names(rna)) {
  rna <- rna %>% mutate(rna_direction = NA_character_)
}

rna <- rna %>%
  mutate(
    legacy_direction_fallback = case_when(
      legacy_rna_regulation == "up" ~ "up",
      legacy_rna_regulation == "down" ~ "down",
      legacy_rna_regulation %in% c("mixed", "unchanged") ~ "no_dominant_direction",
      TRUE ~ "unknown"
    ),
    rna_direction = coalesce(rna_direction, legacy_direction_fallback)
  ) %>%
  select(-legacy_direction_fallback)

say("RNA direction source = ", direction_source)

# ============================================================================
# 3. Build simple baselines on identical modality components
# ============================================================================
say("\n=== STEP 3: calculate baselines ===")

context <- full_join(
  rna,
  chip,
  by = c("tissue", "gene")
) %>%
  mutate(
    rna_present = !is.na(rna_component),
    chip_present = !is.na(chip_component),
    rna_component = coalesce(rna_component, 0),
    chip_component = coalesce(chip_component, 0),
    rna_zero = if_else(rna_present, rna_component, 0),
    chip_zero = if_else(chip_present, chip_component, 0),
    chip_level = case_when(
      !chip_present ~ "none",
      chip_binding_score >= 0.8 ~ "high",
      chip_binding_score >= 0.4 ~ "medium",
      TRUE ~ "low"
    )
  ) %>%
  group_by(tissue) %>%
  mutate(
    rna_rank = percent_rank(rna_zero),
    chip_rank = percent_rank(chip_zero)
  ) %>%
  ungroup() %>%
  mutate(
    # Removes only the biological pattern matrix while retaining
    # the thesis 0.60/0.40 modality balance and 0.70/0.60
    # single-modality penalties.
    Weighted60_40_noPattern = case_when(
      rna_present & chip_present ~ 0.60 * rna_component + 0.40 * chip_component,
      rna_present & !chip_present ~ 0.70 * rna_component,
      !rna_present & chip_present ~ 0.60 * chip_component,
      TRUE ~ 0
    ),

    # Equal-weight late integration; missing modality contributes zero.
    EqualMean_missing0 = 0.50 * rna_zero + 0.50 * chip_zero,

    # Average only available modalities, intentionally removing the
    # incomplete-evidence penalty.
    AvailableMean = case_when(
      rna_present & chip_present ~ 0.50 * (rna_component + chip_component),
      rna_present ~ rna_component,
      chip_present ~ chip_component,
      TRUE ~ 0
    ),

    # Scale-robust rank baselines within tissue.
    MeanRank_missing0 = 0.50 * rna_rank + 0.50 * chip_rank,
    RankProduct_missing0 = sqrt(pmax(rna_rank, 0) * pmax(chip_rank, 0)),

    # Very permissive single-best-evidence baseline.
    MaxEvidence = pmax(rna_zero, chip_zero),

    # Direction-aware concordance proxy based on final modality-level evidence.
    # This is NOT described as raw "DEG ∩ peak"; that would require a different
    # dataset-level binary definition.
    StrictConcordance_proxy = if_else(
      rna_present & chip_present &
        rna_direction == "up" &
        chip_level %in% c("high", "medium"),
      0.60 * rna_component + 0.40 * chip_component,
      0
    )
  )

baseline_cols <- c(
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "AvailableMean",
  "MeanRank_missing0",
  "RankProduct_missing0",
  "MaxEvidence",
  "StrictConcordance_proxy"
)

# Aggregate each baseline exactly like thesis COMPASS:
# max score across non-ESC human-normal tissues.
baseline_gene <- context %>%
  group_by(gene) %>%
  summarise(
    across(all_of(baseline_cols), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Use canonical COMPASS CSV to define the exact universe.
benchmark <- canonical %>%
  left_join(baseline_gene, by = "gene") %>%
  mutate(
    across(all_of(baseline_cols), ~ replace_na(.x, 0))
  )

readr::write_tsv(
  benchmark,
  file.path(OUT_DIR, "03_gene_level_scores_COMPASS_and_baselines.tsv")
)

if (WRITE_FULL_CONTEXT_TABLE) {
  readr::write_csv(
    context,
    file.path(OUT_DIR, "04_context_level_baseline_inputs.csv.gz")
  )
}

coverage <- tibble(
  item = c(
    "canonical_universe",
    "genes_with_any_RNA_baseline_evidence",
    "genes_with_any_ChIP_baseline_evidence",
    "genes_with_both_RNA_and_ChIP_somewhere"
  ),
  n = c(
    nrow(benchmark),
    n_distinct(context$gene[context$rna_present]),
    n_distinct(context$gene[context$chip_present]),
    n_distinct(context$gene[context$rna_present & context$chip_present])
  )
)
readr::write_tsv(coverage, file.path(OUT_DIR, "05_baseline_input_coverage.tsv"))

# ============================================================================
# 4. Benchmark metrics
# ============================================================================
say("\n=== STEP 4: external benchmark ===")

methods <- c("COMPASS", baseline_cols)
metric_rows <- list()
topk_rows <- list()

for (set_name in names(sets)) {
  positives <- sets[[set_name]]
  y <- as.integer(benchmark$gene %in% positives)

  for (m in methods) {
    sc <- benchmark[[m]]

    roc_obj <- pROC::roc(
      y, sc,
      levels = c(0,1),
      direction = "<",
      quiet = TRUE
    )
    auc_v <- as.numeric(pROC::auc(roc_obj))
    ap_v <- average_precision_grouped(y, sc)
    prevalence <- mean(y)

    metric_rows[[length(metric_rows) + 1]] <- tibble(
      reference_set = set_name,
      method = m,
      universe_n = length(y),
      positives_n = sum(y),
      prevalence = prevalence,
      ROC_AUC = auc_v,
      AP = ap_v,
      AP_over_prevalence = ap_v / prevalence
    )

    for (k in TOP_K) {
      z <- topk_metrics(benchmark$gene, y, sc, k)
      topk_rows[[length(topk_rows) + 1]] <- z %>%
        mutate(
          reference_set = set_name,
          method = m,
          .before = 1
        )
    }
  }
}

metrics <- bind_rows(metric_rows)
topks <- bind_rows(topk_rows)

readr::write_tsv(metrics, file.path(OUT_DIR, "06_reference_metrics.tsv"))
readr::write_tsv(topks, file.path(OUT_DIR, "07_topK_metrics.tsv"))

# ============================================================================
# 5. Paired stratified bootstrap: COMPASS minus each baseline
# ============================================================================
say("\n=== STEP 5: paired bootstrap ===")

primary_set <- toupper(PRIMARY_REFERENCE)
positive_genes <- sets[[primary_set]]
y0 <- as.integer(benchmark$gene %in% positive_genes)
pos_idx <- which(y0 == 1)
neg_idx <- which(y0 == 0)

set.seed(SEED)
boot_methods <- baseline_cols
boot <- lapply(
  boot_methods,
  function(x) matrix(
    NA_real_,
    nrow = N_BOOT,
    ncol = 2,
    dimnames = list(NULL, c("delta_auc", "delta_ap"))
  )
)
names(boot) <- boot_methods

for (b in seq_len(N_BOOT)) {
  idx <- c(
    sample(pos_idx, length(pos_idx), replace = TRUE),
    sample(neg_idx, length(neg_idx), replace = TRUE)
  )
  yy <- y0[idx]

  comp_score <- benchmark$COMPASS[idx]
  comp_auc <- as.numeric(
    pROC::auc(
      pROC::roc(
        yy, comp_score,
        levels = c(0,1),
        direction = "<",
        quiet = TRUE
      )
    )
  )
  comp_ap <- average_precision_grouped(yy, comp_score)

  for (m in boot_methods) {
    bb <- benchmark[[m]][idx]
    b_auc <- as.numeric(
      pROC::auc(
        pROC::roc(
          yy, bb,
          levels = c(0,1),
          direction = "<",
          quiet = TRUE
        )
      )
    )
    b_ap <- average_precision_grouped(yy, bb)

    boot[[m]][b, "delta_auc"] <- comp_auc - b_auc
    boot[[m]][b, "delta_ap"] <- comp_ap - b_ap
  }

  if (b %% max(1L, floor(N_BOOT / 10L)) == 0L) {
    message(" bootstrap ", b, "/", N_BOOT)
  }
}

empirical_p <- function(x) {
  left <- (sum(x <= 0, na.rm = TRUE) + 1) / (sum(is.finite(x)) + 1)
  right <- (sum(x >= 0, na.rm = TRUE) + 1) / (sum(is.finite(x)) + 1)
  min(1, 2 * min(left, right))
}

boot_summary <- bind_rows(lapply(boot_methods, function(m) {
  da <- boot[[m]][, "delta_auc"]
  dp <- boot[[m]][, "delta_ap"]

  tibble(
    reference_set = primary_set,
    baseline = m,
    delta_AUC = mean(da, na.rm = TRUE),
    delta_AUC_CI_low = quantile(da, 0.025, na.rm = TRUE),
    delta_AUC_CI_high = quantile(da, 0.975, na.rm = TRUE),
    delta_AUC_p = empirical_p(da),
    delta_AP = mean(dp, na.rm = TRUE),
    delta_AP_CI_low = quantile(dp, 0.025, na.rm = TRUE),
    delta_AP_CI_high = quantile(dp, 0.975, na.rm = TRUE),
    delta_AP_p = empirical_p(dp)
  )
})) %>%
  mutate(
    delta_AUC_FDR = p.adjust(delta_AUC_p, method = "BH"),
    delta_AP_FDR = p.adjust(delta_AP_p, method = "BH")
  )

readr::write_tsv(
  boot_summary,
  file.path(OUT_DIR, "08_bootstrap_COMPASS_vs_baselines.tsv")
)

# ============================================================================
# 6. Method-ranking correlations
# ============================================================================
corr <- suppressWarnings(
  cor(
    benchmark %>% select(all_of(methods)),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
)

corr_df <- as.data.frame(corr) %>%
  tibble::rownames_to_column("method")

readr::write_tsv(
  corr_df,
  file.path(OUT_DIR, "09_method_spearman_correlations.tsv")
)

# ============================================================================
# 7. Compact manuscript-oriented plots
# ============================================================================
p_auc <- metrics %>%
  ggplot(aes(x = method, y = ROC_AUC)) +
  geom_col() +
  facet_wrap(~ reference_set) +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = "ROC-AUC", title = "COMPASS versus simple fusion baselines")

ggsave(
  file.path(OUT_DIR, "figures", "AUC_baselines.pdf"),
  p_auc,
  width = 8.5,
  height = 7
)

p_ap <- metrics %>%
  ggplot(aes(x = method, y = AP_over_prevalence)) +
  geom_col() +
  facet_wrap(~ reference_set) +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Average precision / prevalence",
    title = "Precision-recall lift over background prevalence"
  )

ggsave(
  file.path(OUT_DIR, "figures", "AP_lift_baselines.pdf"),
  p_ap,
  width = 8.5,
  height = 7
)

p_r500 <- topks %>%
  filter(k == PRIMARY_K) %>%
  ggplot(aes(x = method, y = recall)) +
  geom_col() +
  facet_wrap(~ reference_set) +
  coord_flip() +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = paste0("Recall@", PRIMARY_K))

ggsave(
  file.path(OUT_DIR, "figures", "Recall500_baselines.pdf"),
  p_r500,
  width = 8.5,
  height = 7
)

# ============================================================================
# 8. Final report
# ============================================================================
say("\n=== DONE ===")
say("Return the whole folder: ", OUT_DIR)
say("Most important files:")
say("  01_thesis_ROC_reproduction.tsv")
say("  02_sanity_checks.tsv")
say("  03_gene_level_scores_COMPASS_and_baselines.tsv")
say("  05_baseline_input_coverage.tsv")
say("  06_reference_metrics.tsv")
say("  07_topK_metrics.tsv")
say("  08_bootstrap_COMPASS_vs_baselines.tsv")
say("  09_method_spearman_correlations.tsv")

writeLines(log_lines, file.path(OUT_DIR, "00_run_report.txt"))
