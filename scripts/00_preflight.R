#!/usr/bin/env Rscript

# Reproducibility preflight for final journal run.
source("config.R")

required <- c(
  COMPASS_CSV,
  unname(GMT),
  RNA_RDA,
  CHIP_RDA
)

cat("=== COMPASS final journal preflight ===\n")
for (p in required) {
  cat(if (file.exists(p)) "[OK]   " else "[MISS] ", p, "\n", sep="")
}
cat(if (file.exists(RNA_FREQ_RDA)) "[OK]   " else "[WARN] ",
    RNA_FREQ_RDA, "  (optional but recommended)\n", sep="")

if (!all(file.exists(required))) {
  stop("Required inputs are missing.")
}

dir.create("final_reproducibility", showWarnings=FALSE)

# Store hashes so the exact inputs used for the paper are frozen.
md5 <- tools::md5sum(c(required, RNA_FREQ_RDA[file.exists(RNA_FREQ_RDA)]))
write.table(
  data.frame(file=names(md5), md5=unname(md5)),
  "final_reproducibility/input_md5.tsv",
  sep="\t", quote=FALSE, row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  "final_reproducibility/sessionInfo.txt"
)

cat("\nInput hashes and sessionInfo written to final_reproducibility/.\n")
