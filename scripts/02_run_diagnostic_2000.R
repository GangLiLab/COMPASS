#!/usr/bin/env Rscript

# ============================================================================
# COMPASS diagnostic v3
#
# Run AFTER the previous quick benchmark has completed.
#
# Purpose:
# 1) Determine whether the unexpectedly strong AvailableMean / MaxEvidence
#    results are driven by single-modality (especially ChIP-only) maxima.
# 2) Add explicit RNA-only and ChIP-only baselines.
# 3) Re-run the benchmark in a matched "both modalities available" setting,
#    so missing-modality handling cannot dominate the comparison.
# 4) Bootstrap AUC/AP differences for all three Ben-Porath reference sets.
#
# It uses the SAME data/ folder and the existing benchmark_output/ folder.
# No raw sequencing data are required.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
})

DATA_DIR <- "data"
BENCH_DIR <- "benchmark_output"
OUT_DIR <- "diagnostic_output"

GMT <- c(
  PRC2  = file.path(DATA_DIR, "BENPORATH_PRC2_TARGETS.v2026.1.Hs.gmt"),
  SUZ12 = file.path(DATA_DIR, "BENPORATH_SUZ12_TARGETS.v2026.1.Hs.gmt"),
  EED   = file.path(DATA_DIR, "BENPORATH_EED_TARGETS.v2026.1.Hs.gmt")
)

COMPASS_CSV <- file.path(DATA_DIR, "human_normal_all.csv")
CONTEXT_CSV <- file.path(BENCH_DIR, "04_context_level_baseline_inputs.csv.gz")
GENE_TSV <- file.path(BENCH_DIR, "03_gene_level_scores_COMPASS_and_baselines.tsv")

TOP_K <- c(100L, 250L, 500L, 1000L)
N_BOOT <- 2000L
SEED <- 20260727L

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

must_exist <- c(COMPASS_CSV, CONTEXT_CSV, GENE_TSV, unname(GMT))
missing <- must_exist[!file.exists(must_exist)]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"))
}

normalize_gene <- function(x) toupper(trimws(as.character(x)))

normalize_tissue <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- gsub("_", "-", z, fixed = TRUE)
  z <- gsub("[[:space:]]+", "-", z)
  z <- gsub("-+", "-", z)
  z[z == "testis"] <- "testicle"
  z[z == "head-and neck"] <- "head-and-neck"
  z
}

read_gmt_genes <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nchar(ln) > 0][1]
  normalize_gene(strsplit(ln, "\t")[[1]][-(1:2)])
}

auc_rank <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- score[ok]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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
  cp <- cumsum(pos_by)
  cn <- cumsum(n_by)
  precision <- cp / cn
  recall <- cp / npos
  delta_recall <- c(recall[1], diff(recall))
  sum(delta_recall * precision)
}

calc_metrics <- function(df, methods, sets, analysis_name) {
  out <- list()
  top_out <- list()

  for (set_name in names(sets)) {
    y <- as.integer(df$gene %in% sets[[set_name]])
    prevalence <- mean(y)

    for (m in methods) {
      sc <- df[[m]]
      ap <- average_precision_grouped(y, sc)

      out[[length(out) + 1]] <- tibble(
        analysis = analysis_name,
        reference_set = set_name,
        method = m,
        universe_n = nrow(df),
        positives_n = sum(y),
        prevalence = prevalence,
        ROC_AUC = auc_rank(y, sc),
        AP = ap,
        AP_over_prevalence = ap / prevalence
      )

      for (k in TOP_K) {
        kk <- min(k, nrow(df))
        ord <- order(-sc, df$gene)
        idx <- ord[seq_len(kk)]
        tp <- sum(y[idx] == 1)

        top_out[[length(top_out) + 1]] <- tibble(
          analysis = analysis_name,
          reference_set = set_name,
          method = m,
          k = kk,
          true_positives = tp,
          precision = tp / kk,
          recall = tp / sum(y == 1),
          enrichment = (tp / kk) / prevalence
        )
      }
    }
  }

  list(
    metrics = bind_rows(out),
    topk = bind_rows(top_out)
  )
}

bootstrap_deltas <- function(df, methods, sets, analysis_name, B = N_BOOT) {
  baselines <- setdiff(methods, "COMPASS")
  result <- list()
  set.seed(SEED)

  empirical_p <- function(x) {
    p1 <- (sum(x <= 0, na.rm = TRUE) + 1) / (sum(is.finite(x)) + 1)
    p2 <- (sum(x >= 0, na.rm = TRUE) + 1) / (sum(is.finite(x)) + 1)
    min(1, 2 * min(p1, p2))
  }

  for (set_name in names(sets)) {
    y0 <- as.integer(df$gene %in% sets[[set_name]])
    pos <- which(y0 == 1)
    neg <- which(y0 == 0)

    stores <- lapply(
      baselines,
      function(x) matrix(
        NA_real_,
        nrow = B,
        ncol = 2,
        dimnames = list(NULL, c("dAUC", "dAP"))
      )
    )
    names(stores) <- baselines

    for (b in seq_len(B)) {
      idx <- c(
        sample(pos, length(pos), replace = TRUE),
        sample(neg, length(neg), replace = TRUE)
      )
      yy <- y0[idx]
      comp_auc <- auc_rank(yy, df$COMPASS[idx])
      comp_ap <- average_precision_grouped(yy, df$COMPASS[idx])

      for (m in baselines) {
        stores[[m]][b, "dAUC"] <- comp_auc - auc_rank(yy, df[[m]][idx])
        stores[[m]][b, "dAP"] <- comp_ap - average_precision_grouped(yy, df[[m]][idx])
      }
    }

    one_set <- bind_rows(lapply(baselines, function(m) {
      da <- stores[[m]][, "dAUC"]
      dp <- stores[[m]][, "dAP"]

      tibble(
        analysis = analysis_name,
        reference_set = set_name,
        baseline = m,
        delta_AUC = mean(da, na.rm = TRUE),
        delta_AUC_low = quantile(da, 0.025, na.rm = TRUE),
        delta_AUC_high = quantile(da, 0.975, na.rm = TRUE),
        p_AUC = empirical_p(da),
        delta_AP = mean(dp, na.rm = TRUE),
        delta_AP_low = quantile(dp, 0.025, na.rm = TRUE),
        delta_AP_high = quantile(dp, 0.975, na.rm = TRUE),
        p_AP = empirical_p(dp)
      )
    })) %>%
      mutate(
        FDR_AUC = p.adjust(p_AUC, method = "BH"),
        FDR_AP = p.adjust(p_AP, method = "BH")
      )

    result[[length(result) + 1]] <- one_set
  }

  bind_rows(result)
}

# ---------------------------------------------------------------------------
# Load existing outputs
# ---------------------------------------------------------------------------
sets <- lapply(GMT, read_gmt_genes)

gene_scores <- read_tsv(GENE_TSV, show_col_types = FALSE) %>%
  mutate(gene = normalize_gene(gene))

ctx <- read_csv(CONTEXT_CSV, show_col_types = FALSE) %>%
  mutate(
    gene = normalize_gene(gene),
    tissue = normalize_tissue(tissue),
    evidence_type = case_when(
      rna_present & chip_present ~ "both",
      rna_present & !chip_present ~ "rna_only",
      !rna_present & chip_present ~ "chip_only",
      TRUE ~ "none"
    )
  )

# ---------------------------------------------------------------------------
# A. Explicit single-modality baselines
# ---------------------------------------------------------------------------
single_gene <- ctx %>%
  group_by(gene) %>%
  summarise(
    RNAOnly = max(rna_component, na.rm = TRUE),
    ChIPOnly = max(chip_component, na.rm = TRUE),
    .groups = "drop"
  )

full_diag <- gene_scores %>%
  select(
    gene,
    COMPASS,
    Weighted60_40_noPattern,
    EqualMean_missing0,
    AvailableMean,
    MeanRank_missing0,
    RankProduct_missing0,
    MaxEvidence,
    StrictConcordance_proxy
  ) %>%
  left_join(single_gene, by = "gene") %>%
  mutate(
    RNAOnly = replace_na(RNAOnly, 0),
    ChIPOnly = replace_na(ChIPOnly, 0)
  )

full_methods <- c(
  "COMPASS",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "AvailableMean",
  "MaxEvidence",
  "RNAOnly",
  "ChIPOnly",
  "MeanRank_missing0",
  "RankProduct_missing0",
  "StrictConcordance_proxy"
)

full_res <- calc_metrics(full_diag, full_methods, sets, "full_canonical_universe")
write_tsv(full_res$metrics, file.path(OUT_DIR, "01_full_metrics_with_single_modality.tsv"))
write_tsv(full_res$topk, file.path(OUT_DIR, "02_full_topK_with_single_modality.tsv"))

# ---------------------------------------------------------------------------
# B. Provenance of the winning tissue for each baseline
# ---------------------------------------------------------------------------
prov_methods <- c(
  "AvailableMean",
  "MaxEvidence",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "MeanRank_missing0",
  "RankProduct_missing0",
  "StrictConcordance_proxy"
)

prov_rows <- list()

for (m in prov_methods) {
  for (k in TOP_K) {
    top_genes <- gene_scores %>%
      arrange(desc(.data[[m]]), gene) %>%
      slice_head(n = k) %>%
      pull(gene)

    sub <- ctx %>% filter(gene %in% top_genes)
    max_by_gene <- sub %>%
      group_by(gene) %>%
      summarise(max_score = max(.data[[m]], na.rm = TRUE), .groups = "drop")

    winners <- sub %>%
      inner_join(max_by_gene, by = "gene") %>%
      filter(abs(.data[[m]] - max_score) < 1e-12)

    per_gene <- winners %>%
      group_by(gene) %>%
      summarise(
        winner_types = paste(sort(unique(evidence_type)), collapse = "+"),
        .groups = "drop"
      )

    tab <- per_gene %>%
      count(winner_types, name = "n") %>%
      mutate(
        method = m,
        k = k,
        fraction = n / sum(n),
        .before = 1
      )

    prov_rows[[length(prov_rows) + 1]] <- tab
  }
}

provenance <- bind_rows(prov_rows)
write_tsv(provenance, file.path(OUT_DIR, "03_topK_winning_context_provenance.tsv"))

# ---------------------------------------------------------------------------
# C. Matched both-modality context benchmark
# ---------------------------------------------------------------------------
comp_ctx <- fread(COMPASS_CSV)
comp_ctx <- comp_ctx[Gene_biotype == "protein_coding" & Tissue != "esc"]
comp_ctx[, gene := normalize_gene(Gene)]
comp_ctx[, tissue := normalize_tissue(Tissue)]

comp_ctx <- as_tibble(comp_ctx) %>%
  group_by(gene, tissue) %>%
  summarise(COMPASS = max(PTS), .groups = "drop")

matched_ctx <- ctx %>%
  filter(rna_present & chip_present) %>%
  inner_join(comp_ctx, by = c("gene", "tissue"))

matched_methods_context <- c(
  "COMPASS",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "MeanRank_missing0",
  "RankProduct_missing0",
  "MaxEvidence",
  "StrictConcordance_proxy"
)

matched_gene <- matched_ctx %>%
  group_by(gene) %>%
  summarise(
    across(all_of(matched_methods_context), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  )

matched_res <- calc_metrics(
  matched_gene,
  matched_methods_context,
  sets,
  "matched_both_modalities"
)

write_tsv(matched_res$metrics, file.path(OUT_DIR, "04_matched_both_metrics.tsv"))
write_tsv(matched_res$topk, file.path(OUT_DIR, "05_matched_both_topK.tsv"))

matched_summary <- tibble(
  item = c(
    "matched_gene_universe",
    "matched_gene_tissue_rows",
    "PRC2_positives",
    "SUZ12_positives",
    "EED_positives"
  ),
  n = c(
    nrow(matched_gene),
    nrow(matched_ctx),
    sum(matched_gene$gene %in% sets$PRC2),
    sum(matched_gene$gene %in% sets$SUZ12),
    sum(matched_gene$gene %in% sets$EED)
  )
)
write_tsv(matched_summary, file.path(OUT_DIR, "06_matched_both_summary.tsv"))

# ---------------------------------------------------------------------------
# D. Paired bootstrap for all three reference sets
#    Focus on interpretable baselines rather than every possible method.
# ---------------------------------------------------------------------------
full_boot_methods <- c(
  "COMPASS",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "AvailableMean",
  "MaxEvidence",
  "RNAOnly",
  "ChIPOnly",
  "StrictConcordance_proxy"
)

boot_full <- bootstrap_deltas(
  full_diag,
  full_boot_methods,
  sets,
  "full_canonical_universe",
  B = N_BOOT
)

write_tsv(boot_full, file.path(OUT_DIR, "07_bootstrap_all_sets_full.tsv"))

matched_boot_methods <- c(
  "COMPASS",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "MaxEvidence",
  "StrictConcordance_proxy"
)

boot_matched <- bootstrap_deltas(
  matched_gene,
  matched_boot_methods,
  sets,
  "matched_both_modalities",
  B = N_BOOT
)

write_tsv(boot_matched, file.path(OUT_DIR, "08_bootstrap_all_sets_matched_both.tsv"))

# ---------------------------------------------------------------------------
# E. Compact interpretation table
# ---------------------------------------------------------------------------
key <- full_res$metrics %>%
  filter(
    method %in% c(
      "COMPASS",
      "Weighted60_40_noPattern",
      "AvailableMean",
      "MaxEvidence",
      "RNAOnly",
      "ChIPOnly"
    )
  ) %>%
  select(reference_set, method, ROC_AUC, AP, AP_over_prevalence) %>%
  arrange(reference_set, desc(ROC_AUC))

write_tsv(key, file.path(OUT_DIR, "09_key_diagnostic_comparison.tsv"))

cat("\nDiagnostic v3 completed.\n")
cat("Please return the entire diagnostic_output folder.\n")
