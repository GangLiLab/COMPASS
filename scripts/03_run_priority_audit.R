#!/usr/bin/env Rscript

# ============================================================================
# COMPASS priority audit v4
#
# Final focused audit before manuscript freeze.
# No raw sequencing data are reprocessed.
#
# Questions:
# 1) For each method's top-ranked genes, what evidence pattern generated the
#    MAXIMUM score (winning context)?
# 2) Does COMPASS trade some canonical-reference recovery for a higher fraction
#    of convergent RNA-up + H3K27me3-supported winning contexts?
# 3) How different are the top-K candidate sets between methods?
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
})

DATA_DIR <- "data"
BENCH_DIR <- "benchmark_output"
OUT_DIR <- "priority_audit_output"

COMPASS_CSV <- file.path(DATA_DIR, "human_normal_all.csv")
CONTEXT_CSV <- file.path(BENCH_DIR, "04_context_level_baseline_inputs.csv.gz")
GENE_TSV <- file.path(BENCH_DIR, "03_gene_level_scores_COMPASS_and_baselines.tsv")
METRICS_TSV <- file.path(BENCH_DIR, "06_reference_metrics.tsv")

TOP_K <- c(100L, 250L, 500L, 1000L)
EPS <- 1e-10

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

required <- c(COMPASS_CSV, CONTEXT_CSV, GENE_TSV, METRICS_TSV)
if (any(!file.exists(required))) {
  stop("Missing files:\n", paste(required[!file.exists(required)], collapse="\n"))
}

# ------------------------
# Load context information
# ------------------------
ctx <- read_csv(CONTEXT_CSV, show_col_types = FALSE) %>%
  mutate(
    gene = normalize_gene(gene),
    tissue = normalize_tissue(tissue),
    evidence_type = case_when(
      rna_present & chip_present ~ "both",
      rna_present & !chip_present ~ "rna_only",
      !rna_present & chip_present ~ "chip_only",
      TRUE ~ "none"
    ),
    convergent = (
      rna_present & chip_present &
      rna_direction == "up" &
      chip_level %in% c("high", "medium")
    ),
    convergent_highRNA = (
      convergent &
      rna_confidence %in% c("high", "medium")
    )
  )

# COMPASS context score comes from the thesis-canonical CSV.
comp <- fread(COMPASS_CSV)
need <- c("Gene","Gene_biotype","Tissue","PTS")
if (!all(need %in% names(comp))) {
  stop("human_normal_all.csv must contain: ", paste(need, collapse=", "))
}

comp <- comp[
  Gene_biotype == "protein_coding" &
  tolower(Tissue) != "esc"
]
comp[, gene := normalize_gene(Gene)]
comp[, tissue := normalize_tissue(Tissue)]

comp_ctx <- as_tibble(comp) %>%
  group_by(gene, tissue) %>%
  summarise(COMPASS = max(PTS, na.rm=TRUE), .groups="drop")

# Attach evidence features to the exact COMPASS gene × tissue records.
comp_ctx <- comp_ctx %>%
  left_join(
    ctx %>% select(
      gene, tissue, rna_present, chip_present, evidence_type,
      rna_direction, rna_confidence, chip_level,
      convergent, convergent_highRNA
    ) %>% distinct(),
    by=c("gene","tissue")
  )

# Existing method scores
gene_scores <- read_tsv(GENE_TSV, show_col_types=FALSE) %>%
  mutate(gene=normalize_gene(gene))

metrics <- read_tsv(METRICS_TSV, show_col_types=FALSE)

methods <- c(
  "COMPASS",
  "Weighted60_40_noPattern",
  "EqualMean_missing0",
  "AvailableMean",
  "MeanRank_missing0",
  "RankProduct_missing0",
  "MaxEvidence",
  "StrictConcordance_proxy"
)

baseline_methods <- setdiff(methods, "COMPASS")

# ------------------------
# Winning-context function
# ------------------------
get_winners <- function(method) {
  if (method == "COMPASS") {
    dat <- comp_ctx %>%
      transmute(
        gene, tissue,
        score=COMPASS,
        rna_present, chip_present, evidence_type,
        rna_direction, rna_confidence, chip_level,
        convergent, convergent_highRNA
      )
  } else {
    dat <- ctx %>%
      transmute(
        gene, tissue,
        score=.data[[method]],
        rna_present, chip_present, evidence_type,
        rna_direction, rna_confidence, chip_level,
        convergent, convergent_highRNA
      )
  }

  mx <- dat %>%
    group_by(gene) %>%
    summarise(max_score=max(score, na.rm=TRUE), .groups="drop")

  dat %>%
    inner_join(mx, by="gene") %>%
    filter(abs(score-max_score) <= EPS) %>%
    filter(is.finite(score))
}

winner_cache <- lapply(methods, get_winners)
names(winner_cache) <- methods

# ------------------------
# Top-K evidence audit
# ------------------------
summary_rows <- list()
detail_rows <- list()

for (m in methods) {
  win <- winner_cache[[m]]

  for (k in TOP_K) {
    top_genes <- gene_scores %>%
      arrange(desc(.data[[m]]), gene) %>%
      slice_head(n=k) %>%
      pull(gene)

    w <- win %>% filter(gene %in% top_genes)

    # A gene can have tied winning tissues. Aggregate conservatively:
    # "any" indicates whether at least one exact maximum-scoring context
    # satisfies the evidence property.
    pg <- w %>%
      group_by(gene) %>%
      summarise(
        winning_context_n=n(),
        winner_any_both=any(evidence_type=="both", na.rm=TRUE),
        winner_any_rna_only=any(evidence_type=="rna_only", na.rm=TRUE),
        winner_any_chip_only=any(evidence_type=="chip_only", na.rm=TRUE),
        winner_any_convergent=any(convergent %in% TRUE, na.rm=TRUE),
        winner_any_convergent_highRNA=any(convergent_highRNA %in% TRUE, na.rm=TRUE),
        winning_tissues=paste(sort(unique(tissue)), collapse=";"),
        .groups="drop"
      )

    # Preserve genes whose winning context metadata failed to join.
    pg <- tibble(gene=top_genes) %>%
      left_join(pg, by="gene") %>%
      mutate(
        across(
          c(winner_any_both, winner_any_rna_only, winner_any_chip_only,
            winner_any_convergent, winner_any_convergent_highRNA),
          ~replace_na(.x, FALSE)
        ),
        winning_context_n=replace_na(winning_context_n,0L),
        winning_tissues=replace_na(winning_tissues,"")
      )

    summary_rows[[length(summary_rows)+1]] <- tibble(
      method=m,
      k=k,
      top_genes=nrow(pg),
      winner_both_fraction=mean(pg$winner_any_both),
      winner_single_only_fraction=mean(!pg$winner_any_both),
      winner_convergent_fraction=mean(pg$winner_any_convergent),
      winner_convergent_highRNA_fraction=mean(pg$winner_any_convergent_highRNA),
      metadata_join_failure_fraction=mean(pg$winning_context_n==0)
    )

    if (k == 500L) {
      detail_rows[[length(detail_rows)+1]] <- pg %>%
        mutate(method=m, .before=1)
    }
  }
}

priority_summary <- bind_rows(summary_rows)
priority_detail <- bind_rows(detail_rows)

write_tsv(priority_summary, file.path(OUT_DIR,"01_topK_winning_context_quality.tsv"))
write_tsv(priority_detail, file.path(OUT_DIR,"02_top500_winning_context_detail.tsv"))

# ------------------------
# "Any concordant context" (less strict than winning-context concordance)
# ------------------------
any_concordance <- ctx %>%
  group_by(gene) %>%
  summarise(
    has_any_both=any(evidence_type=="both"),
    has_any_convergent=any(convergent),
    has_any_convergent_highRNA=any(convergent_highRNA),
    .groups="drop"
  )

any_rows <- list()
for (m in methods) {
  for (k in TOP_K) {
    tg <- gene_scores %>%
      arrange(desc(.data[[m]]), gene) %>%
      slice_head(n=k) %>%
      select(gene) %>%
      left_join(any_concordance, by="gene") %>%
      mutate(
        across(starts_with("has_"), ~replace_na(.x,FALSE))
      )

    any_rows[[length(any_rows)+1]] <- tibble(
      method=m, k=k,
      any_both_fraction=mean(tg$has_any_both),
      any_convergent_fraction=mean(tg$has_any_convergent),
      any_convergent_highRNA_fraction=mean(tg$has_any_convergent_highRNA)
    )
  }
}
write_tsv(bind_rows(any_rows), file.path(OUT_DIR,"03_topK_any_concordant_context.tsv"))

# ------------------------
# Top-K Jaccard overlap with COMPASS
# ------------------------
jrows <- list()
for (k in TOP_K) {
  comp_set <- gene_scores %>%
    arrange(desc(COMPASS), gene) %>%
    slice_head(n=k) %>%
    pull(gene) %>%
    unique()

  for (m in baseline_methods) {
    s <- gene_scores %>%
      arrange(desc(.data[[m]]), gene) %>%
      slice_head(n=k) %>%
      pull(gene) %>%
      unique()

    inter <- length(intersect(comp_set,s))
    union <- length(union(comp_set,s))

    jrows[[length(jrows)+1]] <- tibble(
      method=m,
      k=k,
      overlap_n=inter,
      jaccard=inter/union
    )
  }
}
write_tsv(bind_rows(jrows), file.path(OUT_DIR,"04_topK_jaccard_vs_COMPASS.tsv"))

# ------------------------
# Performance / evidence-quality trade-off table
# ------------------------
q500 <- priority_summary %>%
  filter(k==500) %>%
  select(
    method,
    winner_both_fraction,
    winner_single_only_fraction,
    winner_convergent_fraction,
    winner_convergent_highRNA_fraction
  )

tradeoff <- metrics %>%
  inner_join(q500, by="method") %>%
  arrange(reference_set, desc(ROC_AUC))

write_tsv(tradeoff, file.path(OUT_DIR,"05_performance_vs_evidence_quality.tsv"))

cat("\nPriority audit completed.\n")
cat("Return the entire priority_audit_output folder.\n")
